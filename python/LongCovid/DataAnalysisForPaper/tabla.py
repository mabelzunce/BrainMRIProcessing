import pandas as pd 
import numpy as np 
from scipy.stats import chi2_contingency, kruskal, ttest_ind



def obtain_values_tests(control_values, covid_values, test, test_puntaje_bruto, test_valor_referencia, test_rendimiento, output):
    
    output.write(f'{test}\n')

    # COVID 
    puntaje_bruto_mean_covid = covid_values[test_puntaje_bruto].mean()
    puntaje_bruto_std_covid = covid_values[test_puntaje_bruto].std()

    percentil_mean_covid = covid_values[test_valor_referencia].mean()
    percentil_std_covid = covid_values[test_valor_referencia].std()

    total_covid = len(covid_values)
    normales_total_covid = len(covid_values.loc[covid_values[test_rendimiento] == 'normal'])
    porcentaje_normales_covid = normales_total_covid/total_covid * 100

    levemente_disminuidos_total_covid = len(covid_values.loc[covid_values[test_rendimiento] == 'levemente disminuido'])
    porcentaje_lev_disminuido_covid = levemente_disminuidos_total_covid/total_covid * 100
    
    disminuidos_total_covid = len(covid_values.loc[covid_values[test_rendimiento] == 'disminuido'])
    porcentaje_disminuido_covid = disminuidos_total_covid/total_covid * 100


    

    output.write(f'COVID Puntaje bruto:{puntaje_bruto_mean_covid:.1f}+{puntaje_bruto_std_covid:.1f}\n')
    output.write(f'COVID percentil={percentil_mean_covid:.1f}+{percentil_std_covid:.1f}\n')
    output.write(f'COVID normales:{normales_total_covid} ({porcentaje_normales_covid:.1f}%),levemente disminuidos={levemente_disminuidos_total_covid}({porcentaje_lev_disminuido_covid:.1f}%), disminuidos={disminuidos_total_covid}({porcentaje_disminuido_covid:.1f}%)\n')
    

    # Control 
    puntaje_bruto_mean_control = control_values[test_puntaje_bruto].mean()
    puntaje_bruto_std_control = control_values[test_puntaje_bruto].std()

    percentil_mean_control =control_values[test_valor_referencia].mean()
    percentil_std_control = control_values[test_valor_referencia].std()

    total_control= len(control_values)
    normales_total_control = len(control_values.loc[control_values[test_rendimiento] == 'normal'])
    porcentaje_normales_control = normales_total_control/total_control * 100

    levemente_disminuidos_total_control = len(control_values.loc[control_values[test_rendimiento] == 'levemente disminuido'])
    porcentaje_lev_disminuido_control = levemente_disminuidos_total_control/total_control * 100

    disminuidos_total_control = len(control_values.loc[control_values[test_rendimiento] == 'disminuido'])
    porcentaje_disminuido_control = disminuidos_total_control/total_control * 100

    output.write(f'Control Puntaje bruto:{puntaje_bruto_mean_control:.1f}+{puntaje_bruto_std_control:.1f}\n')
    output.write(f'Control percentil={percentil_mean_control:.1f}+{percentil_std_control:.1f}\n')
    output.write(f'Control normales:{normales_total_control} ({porcentaje_normales_control:.1f}%), levemente disminuidos={levemente_disminuidos_total_control}({porcentaje_lev_disminuido_control:.1f}%), disminuidos={disminuidos_total_control}({porcentaje_disminuido_control:.1f}%)\n')

    # Estadística
    # Percentil Kruskal Wallis


    # Chi square (Categorias: normales/lev_disminuido/disminuido)
    # Tabla de contigencia 
    covid_counts = covid_values[test_rendimiento].value_counts()
    control_counts = control_values[test_rendimiento].value_counts()

    contingency_table = pd.DataFrame([covid_counts, control_counts], index=['COVID', 'Control']).fillna(0).astype(int)
    
    # Realizar el test Chi-Square
    chi2, p_chi_square, dof, expected = chi2_contingency(contingency_table)
    
    
    # Kruskal Wallis
    # Remove NaN values
    control_data_clean = control_values[test_valor_referencia][~np.isnan(control_values[test_valor_referencia])]
    covid_data_clean = covid_values[test_valor_referencia][~np.isnan(covid_values[test_valor_referencia])]
    stat, p_kruskal = kruskal(control_data_clean, covid_data_clean)
   
    output.write(f'P value Chi Square: {p_chi_square:.2f}\n')
    output.write(f'P value Kruskal Wllis: {p_kruskal:.2f}\n')
    output.write(f'\n')

def obtain_values_questionnaires(control_values, covid_values, test, output):
    
    categorias = covid_values[test].unique()


    total_controles = len(control_values)
    total_covid = len(covid_values)

    output.write(f'{test}\n')
    for categoria in categorias:
        total_categoria_covid = len(covid_values.loc[covid_values[test] == categoria])
        porcentaje_categoria_covid = total_categoria_covid / total_covid * 100
        output.write(f'COVID {categoria}: {total_categoria_covid}({porcentaje_categoria_covid:.1f})%\n')
        
        total_categoria_control = len(control_values.loc[control_values[test] == categoria])
        porcentaje_categoría_control = total_categoria_control / total_controles * 100
        output.write(f'Control {categoria}: {total_categoria_control}({porcentaje_categoría_control:.1f})%\n')
        output.write('\n')

    # Chi square (Categorias: normales/lev_disminuido/disminuido)
    # Tabla de contigencia
    covid_counts = covid_values[test].value_counts()
    control_counts = control_values[test].value_counts()

    contingency_table = pd.DataFrame([covid_counts, control_counts], index=['COVID', 'Control']).fillna(0).astype(
        int)

    # Realizar el test Chi-Square
    chi2, p_chi_square, dof, expected = chi2_contingency(contingency_table)

    output.write(f'P chi square:{p_chi_square}\n')
    output.write('\n')
    return 

path_respuestas_cuestionarios = '/mnt/d87cc26d-5470-443c-81c1-e09b68ee4730/Descargas/RespuestasCuestionarioEvaluacionCognitivaResonancia_2.xlsx'
path_output_test_txt = '/mnt/d87cc26d-5470-443c-81c1-e09b68ee4730/Descargas/output_tests.txt'
path_output_questionnaires = '/mnt/d87cc26d-5470-443c-81c1-e09b68ee4730/Descargas/output_questionnaires.txt'

exclude_subject_cognitive = ['CP0011']
exclude_subject_general = ['CP0015']

df_respuestas_cuestionarios = pd.read_excel(path_respuestas_cuestionarios)

# Filter out the subjects
df_respuestas_cuestionarios_filtered = df_respuestas_cuestionarios[~df_respuestas_cuestionarios['ID'].isin(exclude_subject_general)]
df_filtered_cognitive_test = df_respuestas_cuestionarios_filtered[~df_respuestas_cuestionarios_filtered['ID'].isin(exclude_subject_cognitive)]

COVID_values = df_filtered_cognitive_test.loc[df_filtered_cognitive_test['Grupo'] == "COVID"]
control_values = df_filtered_cognitive_test.loc[df_filtered_cognitive_test['Grupo'] == "CONTROL"]

output = open(path_output_test_txt,'w')
output_questionnaires = open(path_output_questionnaires,'w')

obtain_values_tests(control_values, COVID_values, 'TMT-A','TMTA-puntaje bruto', 
                    'TMTA-valor de referencia', 'TMTA-rendimiento', output)

obtain_values_tests(control_values, COVID_values, 'TMT-B','TMTB-puntaje bruto', 
                    'TMTB-valor de referencia', 'TMTB-rendimiento', output)

obtain_values_tests(control_values, COVID_values, 'WMS-R-directos',
                    'WMS-R-directos puntaje', 'WMS-R-directos valor de referencia', 
                    'WMS-R-directos rendimiento', output)

obtain_values_tests(control_values, COVID_values, 'WMS-R-inversos',
                    'WMS-R-inversos puntaje', 'WMS-R-inversos valor de referencia', 
                    'WMS-R-inversos rendimiento', output)

obtain_values_tests(control_values, COVID_values, 'STROOP-P',
                    'STROOP-P puntaje', 'STROOP-P valor de referencia', 
                    'STROOP-P rendimiento', output)

obtain_values_tests(control_values, COVID_values, 'STROOP-C',
                    'STROOP-C puntaje', 'STROOP-C valor de referencia', 
                    'STROOP-C rendimiento', output)

obtain_values_tests(control_values, COVID_values, 'STROOP-P/C',
                    'STROOP-P/C puntaje bruto', 'STROOP-P/C valor de referencia', 
                    'STROOP-P/C rendimiento', output)

obtain_values_tests(control_values, COVID_values, 'STROOP-interf',
                    'STROOP-interf puntaje', 'STROOP-interf valor de referencia', 
                    'STROOP-interf rendimiento', output)

obtain_values_tests(control_values, COVID_values, 'MOCA',
                    'MOCA puntaje bruto', 'MOCA valor de referencia', 
                    'MOCA rendimiento', output)
# questionnaires
COVID_values_questionnaires = df_respuestas_cuestionarios_filtered.loc[df_respuestas_cuestionarios_filtered['Grupo'] == "COVID"]
control_values_questionnaires = df_respuestas_cuestionarios_filtered.loc[df_respuestas_cuestionarios_filtered['Grupo'] == "CONTROL"]

obtain_values_questionnaires(control_values_questionnaires,  COVID_values_questionnaires, 'CATEGORÍA FAS', output_questionnaires)
obtain_values_questionnaires(control_values_questionnaires,  COVID_values_questionnaires, 'CATEGORÍA IPAQ', output_questionnaires)
obtain_values_questionnaires(control_values_questionnaires,  COVID_values_questionnaires, 'PSQI-categoria', output_questionnaires)
obtain_values_questionnaires(control_values_questionnaires,  COVID_values_questionnaires, 'EQ-5D-5L-Mobility', output_questionnaires)
obtain_values_questionnaires(control_values_questionnaires,  COVID_values_questionnaires, 'EQ-5D-5L-Self-care ', output_questionnaires)
obtain_values_questionnaires(control_values_questionnaires,  COVID_values_questionnaires, 'EQ-5D-5L-Activity', output_questionnaires)
obtain_values_questionnaires(control_values_questionnaires,  COVID_values_questionnaires, 'EQ-5D-5L-Pain', output_questionnaires)
obtain_values_questionnaires(control_values_questionnaires,  COVID_values_questionnaires, 'EQ-5D-5L-Anxiety', output_questionnaires)

# Eq Vas
eq_vas_value_covid_mean = COVID_values_questionnaires['EQ_VAS'].mean()
eq_vas_value_covid_std = COVID_values_questionnaires['EQ_VAS'].std()

eq_vas_value_control_mean = control_values_questionnaires['EQ_VAS'].mean()
eq_vas_value_control_std = control_values_questionnaires['EQ_VAS'].std()

output_questionnaires.write(f'Control EQ_VAS:{eq_vas_value_control_mean:.1f}+{eq_vas_value_control_std:.1f}\n')
output_questionnaires.write(f'COVID EQ_VAS:{eq_vas_value_covid_mean:.1f}+{eq_vas_value_covid_std:.1f}\n')

control_data_eq_vas_clean = control_values_questionnaires['EQ_VAS'][~np.isnan(control_values_questionnaires['EQ_VAS'])]
covid_data_eq_vas_clean = COVID_values_questionnaires['EQ_VAS'][~np.isnan(COVID_values_questionnaires['EQ_VAS'])]

stat, p_kruskal = kruskal(control_data_eq_vas_clean, covid_data_eq_vas_clean)

output_questionnaires.write(f'P value Kruskal Wllis: {p_kruskal}\n')