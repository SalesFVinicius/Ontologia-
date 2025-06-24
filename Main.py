# -*- coding: utf-8 -*-
"""
Created on Sun Jun  1 09:11:21 2025

@author: VINICIUSSALES
"""


from Spectrum_Deep import Picos_Abs
from Ontologia import Ontologia


Function = Picos_Abs()
data = Function.return_p()
onto = Ontologia()

for dat in data:
    
    nome_amostra = dat['nome_amostra']
    bandas_detectadas = dat['bandas_detectadas']
    nome_afloramento = dat['nome_afloramento']
    onto.registro(nome_amostra, nome_afloramento, bandas_detectadas,threshold=4)

onto.salvar_ontologia()


