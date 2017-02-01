import pandas as p

from os import listdir; 

output_dir = '../january-replication/by-draw'

output_list = [f for f in listdir(output_dir) 
               if f.startswith('output')]

import json


filtered_states = p.DataFrame()
smoothed_states = p.DataFrame()
smoothed_shocks = p.DataFrame()
filtered_shocks = p.DataFrame()
j = 0
from tqdm import tqdm 
for output in tqdm(output_list):
    
    data = json.load(open(output_dir+'/'+output))

    for sim in data['output']:

        smooth_i = p.DataFrame(data['output'][sim]['smoothed_states'])
        filter_i = p.DataFrame(data['output'][sim]['filtered_states'])

        smooth_i['trial'] = j
        filter_i['trial'] = j

        smoothed_states = smoothed_states.append(smooth_i)
        filtered_states = filtered_states.append(filter_i)

        smooth_i = p.DataFrame(data['output'][sim]['smoothed_shocks'])
        filter_i = p.DataFrame(data['output'][sim]['filtered_shocks'])

        smooth_i['trial'] = j
        filter_i['trial'] = j

        smoothed_shocks = smoothed_shocks.append(smooth_i)
        filtered_shocks = filtered_shocks.append(filter_i)



        j = j+1

#print(output_list)