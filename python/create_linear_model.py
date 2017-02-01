from dsge.DSGE import DSGE
from dsge.translate import translate
#model_file = '/msu/scratch3/m1cjg01/aer_revision_ed/final_code/input/model_v5_2'
base_model_file = '/msu/scratch3/m1cjg01/aer_revision_ed/python/linearized_model_5obs_gamxhp_JPT_priors.yaml'
model = DSGE.read(model_file)


measurement_error = [0, 0.2, 0.5]
nobs = 5, 7

output_dir = '/msu/scratch3/m1cjg01/aer_revision_ed/final_code/generated_linear_model/'
translate(model, output_dir=output_dir)

with open(output_dir+'/model/trans.txt') as f:
    lines = f.readlines()

# these normal distributions are actually truncated
lines[3] = '1, 0, 999, 1\n'
lines[6] = '1, 0, 999999, 1\n'
lines[7] = '1, 0, 999999, 1\n'
lines[15] = '1, 0, 999, 1\n'
print (''.join(lines))
with open(output_dir+'/model/trans.txt', 'w') as f:
    f.write(''.join(lines))
