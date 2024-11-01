import pandas as pd
import numpy as np
from scipy.stats import chi2_contingency, kruskal, ttest_ind

# Función para aplicar la condición
def education(nivel_educativo):
    if nivel_educativo == 'Posgrado':
        return 7
    elif nivel_educativo == 'Terciario completo':
        return 5
    elif nivel_educativo == 'Primario completo':
        return 1
    elif nivel_educativo == 'Universitario completo':
        return 6
    elif nivel_educativo == 'Secundario completo':
        return 3
    elif nivel_educativo == 'Universitario incompleto':
        return 4
    elif nivel_educativo == 'Terciario incompleto':
        return 4
    elif nivel_educativo == "Secundario incompleto":
        return 2
    elif nivel_educativo == "Primario incompleto":
        return 0
    else:
        return np.nan

path_respuestas_cuestionarios = '/mnt/d87cc26d-5470-443c-81c1-e09b68ee4730/COVID/RespuestasCuestionarioEvaluacionCognitivaResonancia_2.xlsx'
path_respuestas_cuestionarios_out ='/mnt/d87cc26d-5470-443c-81c1-e09b68ee4730/COVID/RespuestasCuestionarioEvaluacionCognitivaResonancia_3.xlsx'

df_respuestas_cuestionarios = pd.read_excel(path_respuestas_cuestionarios)
df_respuestas_cuestionarios['Education_level_cine_11'] = df_respuestas_cuestionarios["Nivel de estudios"].apply(education)

df_respuestas_cuestionarios.to_excel(path_respuestas_cuestionarios_out, index=False)
