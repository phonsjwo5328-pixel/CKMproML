# -*- coding: utf-8 -*-
"""
Step 03: Build the FINAL external-validation dataset.
Design (per published NHANES-CKM ML literature, e.g. Redox Biol 2025;39695774;
Sci Total Environ 2025;40609416; Cardiorenal Med 2025;40690909):
  Population : NHANES 2017-March 2020 adults >=20, non-pregnant, no cancer
  Outcome    : advanced CKM syndrome (stage 3/4) at survey - cross-sectional
               surrogate for "CKM progression", standard practice in NHANES
               CKM studies. 1 = CKM stage 3/4; 0 = stage 0-2.
  Predictors : the 12 CHARLS-model variables in SI units, model naming.
               bodily_pain mapped by composite pain definition:
                 1 = arthritis dx (MCQ160A==1) OR abdominal pain past 12m
                     (MCQ520==1) OR ever seen doctor about pain (MCQ540==1)
  Units      : creatinine/uric acid µmol/L, hemoglobin g/L, HbA1c %,
               HDL/TG/glucose mmol/L, hs-CRP mg/L, SBP mmHg, BMI kg/m²
Outputs:
  wbyz/output/external_validation_final.csv      (analysis-ready)
  wbyz/output/external_validation_codebook.csv   (variable mapping)
  wbyz/output/extraction_log_YYYYMMDD.txt       (this run's numbers)
"""
import sys, os, warnings, datetime
import numpy as np
import pandas as pd

sys.stdout.reconfigure(encoding='utf-8')
warnings.filterwarnings('ignore')

RAW = os.environ.get('NHANES_RAW', 'raw')   # folder holding CDC P-cycle XPT files
OUT = os.environ.get('NHANES_OUT', 'data_public')
os.makedirs(OUT, exist_ok=True)
log = []

def P(msg):
    print(msg)
    log.append(str(msg))

def rd(name, cols):
    df = pd.read_sas(os.path.join(RAW, name + '.xpt'), format='xport')
    keep = [c for c in cols if c in df.columns]
    return df[keep].copy()

P(f'=== NHANES external validation build @ {datetime.date.today()} ===')

# ---------- load ----------
demo = rd('P_DEMO', ['SEQN','RIAGENDR','RIDAGEYR','RIDRETH3','DMDEDUC2','DMDMARTZ',
                     'RIDEXPRG','WTINTPRP','WTMECPRP','SDMVPSU','SDMVSTRA'])
bmx  = rd('P_BMX', ['SEQN','BMXBMI','BMXWAIST'])
bpx  = rd('P_BPXO', ['SEQN','BPXOSY1','BPXOSY2','BPXOSY3','BPXODI1','BPXODI2','BPXODI3'])
ghb  = rd('P_GHB', ['SEQN','LBXGH'])
glu  = rd('P_GLU', ['SEQN','LBDGLUSI'])
hdl  = rd('P_HDL', ['SEQN','LBDHDDSI'])
trig = rd('P_TRIGLY', ['SEQN','WTSAFPRP','LBDTRSI'])
bio  = rd('P_BIOPRO', ['SEQN','LBXSCR','LBXSUA'])
crp  = rd('P_HSCRP', ['SEQN','LBXHSCRP'])
alb  = rd('P_ALB_CR', ['SEQN','URDACT'])
smq  = rd('P_SMQ', ['SEQN','SMQ020','SMQ040'])
alq  = rd('P_ALQ', ['SEQN','ALQ111','ALQ121'])
diq  = rd('P_DIQ', ['SEQN','DIQ010'])
mcq  = rd('P_MCQ', ['SEQN','MCQ160B','MCQ160C','MCQ160D','MCQ160E','MCQ160M','MCQ220',
                    'MCQ160A','MCQ520','MCQ540'])

df = demo
for x in [bmx, bpx, ghb, glu, hdl, trig, bio, crp, alb, smq, alq, diq, mcq]:
    df = df.merge(x, on='SEQN', how='left')
P(f'merged raw: {df.shape[0]} participants, {df.shape[1]} variables')

# ---------- exclusions ----------
n_all = len(df)
df = df[df['RIDAGEYR'] >= 20];  n_20 = len(df)
df = df[df['RIDEXPRG'] != 1];  n_preg = len(df)
df = df[~df['MCQ220'].isin([1, 7, 9])];  n_cancer = len(df)
P(f'exclusions: all {n_all} -> 20+ {n_20} -> non-pregnant {n_preg} -> no-cancer {n_cancer}')

# ---------- BP ----------
sbp_c = ['BPXOSY1','BPXOSY2','BPXOSY3']
dbp_c = ['BPXODI1','BPXODI2','BPXODI3']
df['SBP'] = df[sbp_c].where(df[sbp_c] > 0).mean(axis=1)
df['DBP'] = df[dbp_c].where(df[dbp_c] > 0).mean(axis=1)

# ---------- eGFR (2021 CKD-EPI, no race) ----------
scr = df['LBXSCR']
female = (df['RIAGENDR'] == 2)
egfr = pd.Series(np.nan, index=df.index)
for cond, kappa, alpha, sf in [(female, 0.7, -0.241, 1.012), (~female, 0.9, -0.302, 1.0)]:
    s = scr[cond]
    egfr[cond] = 142 * np.minimum(s / kappa, 1) ** alpha * \
                 np.maximum(s / kappa, 1) ** (-1.200) * (0.9938 ** df['RIDAGEYR'][cond]) * sf
df['e_gfr'] = egfr

# ---------- behaviors ----------
df['smoking_status'] = np.where(df['SMQ020'] == 1, np.where(df['SMQ040'].isin([1, 2]), 1, 0),
                        np.where(df['SMQ020'] == 2, 0, np.nan))
df['alcohol_consumption'] = np.where(df['ALQ111'] == 1, np.where(df['ALQ121'].between(1, 10), 1, 0),
                             np.where(df['ALQ111'] == 2, 0, np.nan))

# ---------- composite bodily pain (user-confirmed mapping) ----------
df['bodily_pain'] = np.where(
    (df['MCQ160A'] == 1) | (df['MCQ520'] == 1) | (df['MCQ540'] == 1), 1, 0)

# ---------- CKM staging (AHA 2023 US) ----------
cvd = (df['MCQ160B'].fillna(2) == 1) | (df['MCQ160C'].fillna(2) == 1) | \
      (df['MCQ160D'].fillna(2) == 1) | (df['MCQ160E'].fillna(2) == 1) | \
      (df['MCQ160M'].fillna(2) == 1)
adv_ckd = (df['e_gfr'] < 30)
albu = (df['URDACT'] >= 30)
stage34 = cvd | adv_ckd | albu

dm = (df['LBDGLUSI'] >= 7.0) | (df['LBXGH'] >= 6.5) | (df['DIQ010'].fillna(2) == 1)
ckd23 = (df['e_gfr'] >= 30) & (df['e_gfr'] < 60)
tg_high = (df['LBDTRSI'] >= 1.7)
stage2 = dm | ckd23 | tg_high

over = (df['BMXBMI'] >= 25) | ((df['RIAGENDR'] == 1) & (df['BMXWAIST'] >= 102)) | \
       ((df['RIAGENDR'] == 2) & (df['BMXWAIST'] >= 88))
predm = df['LBDGLUSI'].between(5.6, 6.9) | df['LBXGH'].between(5.7, 6.4)
dys = (df['LBDHDDSI'] < np.where(female, 1.3, 1.0)) | (df['LBDTRSI'] >= 1.7)
prehy = (df['SBP'] >= 130) | (df['DBP'] >= 80)
stage1 = over | predm | dys | prehy | (df['smoking_status'] == 1)

stage = np.full(len(df), np.nan)
stage[stage34.values] = 3
stage[(~stage34) & stage2.values] = 2
stage[(~stage34) & (~stage2) & stage1.values] = 1
stage[(~stage34) & (~stage2) & (~stage1)] = 0
df['ckm_stage'] = stage
P('\nCKM staging distribution (n=%d):' % len(df))
P(df['ckm_stage'].value_counts(dropna=False).sort_index().to_string())

# ---------- outcome: advanced CKM (stage 3/4) ----------
df['ckm_y'] = np.where(df['ckm_stage'] == 3, 1, np.where(df['ckm_stage'].isna(), np.nan, 0))

# ---------- final 12-variable frame (model naming, SI units) ----------
out = pd.DataFrame({
    'SEQN': df['SEQN'],
    'ckm_y': df['ckm_y'],
    'age': df['RIDAGEYR'],
    'alcohol_consumption': df['alcohol_consumption'],
    'bodily_pain': df['bodily_pain'],
    'body_mass_index': df['BMXBMI'],
    'c_reactive_protein': df['LBXHSCRP'],
    'fasting_glucose': df['LBDGLUSI'],
    'hb_a1c': df['LBXGH'],
    'hdl_c': df['LBDHDDSI'],
    'serum_creatinine_umol_l': df['LBXSCR'] * 88.4,
    'smoking_status': df['smoking_status'],
    'systolic_blood_pressure': df['SBP'],
    'triglyceride': df['LBDTRSI'],
    # ---- transparency / sensitivity-analysis fields (NOT model inputs) ----
    'ckm_stage': df['ckm_stage'],
    'e_gfr': df['e_gfr'],
    'gender': df['RIAGENDR'],
    'race_ethnicity': df['RIDRETH3'],
    'fasting_weight': df['WTSAFPRP'] if 'WTSAFPRP' in df.columns else np.nan,
    'mec_weight': df['WTMECPRP'],
    'sdmvpsu': df['SDMVPSU'],
    'sdmvstra': df['SDMVSTRA'],
})

# complete-case on the 12 predictors + outcome
pred12 = ['age','alcohol_consumption','bodily_pain','body_mass_index','c_reactive_protein',
          'fasting_glucose','hb_a1c','hdl_c','serum_creatinine_umol_l','smoking_status',
          'systolic_blood_pressure','triglyceride']
complete = out.dropna(subset=pred12 + ['ckm_y'])
P(f'\ncomplete-case (12 predictors + outcome): {len(complete)} of {len(out)}')
P(f'outcome: stage 3/4 = {(complete["ckm_y"]==1).sum()}, stage 0-2 = {(complete["ckm_y"]==0).sum()}')
P(f'prevalence of advanced CKM: {(complete["ckm_y"]==1).mean()*100:.1f}%')

# ---------- save ----------
complete.to_csv(os.path.join(OUT, 'external_validation_final.csv'),
                index=False, encoding='ascii', errors='replace')
P('\nsaved: external_validation_final.csv ' + str(complete.shape))
out.to_csv(os.path.join(OUT, 'external_validation_with_missing.csv'),
           index=False, encoding='ascii', errors='replace')
P('saved: external_validation_with_missing.csv ' + str(out.shape))

# ---------- codebook ----------
codebook = [
    ('SEQN','NHANES participant ID','-',''),
    ('ckm_y','Outcome: advanced CKM (stage 3/4), cross-sectional surrogate','0/1','AHA 2023 staging at survey'),
    ('age','Age','years','RIDAGEYR'),
    ('alcohol_consumption','Current drinker','0/1','ALQ111==1 & ALQ121 in 1..10'),
    ('bodily_pain','Composite pain (arthritis dx / abdominal pain 12m / seen doctor for pain)','0/1','MCQ160A | MCQ520 | MCQ540'),
    ('body_mass_index','BMI','kg/m2','BMXBMI'),
    ('c_reactive_protein','hs-CRP','mg/L','LBXHSCRP'),
    ('fasting_glucose','Fasting glucose (fasting subsample)','mmol/L','LBDGLUSI'),
    ('hb_a1c','HbA1c','%','LBXGH'),
    ('hdl_c','HDL cholesterol','mmol/L','LBDHDDSI'),
    ('serum_creatinine_umol_l','Serum creatinine','umol/L','LBXSCR x 88.4'),
    ('smoking_status','Current smoker','0/1','SMQ020==1 & SMQ040 in {1,2}'),
    ('systolic_blood_pressure','SBP mean of up to 3 oscillometric readings','mmHg','BPXOSY1-3'),
    ('triglyceride','Fasting triglycerides (fasting subsample)','mmol/L','LBDTRSI'),
    ('ckm_stage','CKM stage 0-3','0/1/2/3','AHA 2023 US operationalization'),
    ('e_gfr','eGFR 2021 CKD-EPI creatinine no race','mL/min/1.73m2','computed from LBXSCR'),
    ('gender','Sex','1=M,2=F','RIAGENDR'),
    ('race_ethnicity','Race/ethnicity','1..7','RIDRETH3'),
    ('fasting_weight','Fasting subsample 2-yr weight','-','WTSAFPRP'),
    ('mec_weight','MEC exam 2-yr weight','-','WTMECPRP'),
    ('sdmvpsu','PSU (variance estimation)','-','SDMVPSU'),
    ('sdmvstra','Stratum (variance estimation)','-','SDMVSTRA'),
]
cb = pd.DataFrame(codebook, columns=['variable','description','unit','NHANES_source_or_rule'])
cb.to_csv(os.path.join(OUT, 'external_validation_codebook.csv'), index=False, encoding='ascii', errors='replace')
P('saved: external_validation_codebook.csv')

with open(os.path.join(OUT, 'extraction_log_%s.txt' % datetime.date.today().strftime('%Y%m%d')), 'w', encoding='utf-8') as f:
    f.write('\n'.join(log))
P('saved: extraction log')