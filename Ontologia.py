from owlready2 import get_ontology, Thing, ObjectProperty, types
import tkinter as tk
from tkinter import filedialog
import numpy as np

class Ontologia:
    def __init__(self):
        root = tk.Tk()
        root.withdraw()
        print('Import Ontologia')                  
        self.ontology_path = filedialog.askopenfilename(
            parent=None, 
            title="Select file", 
            filetypes=[("rdf", '*.rdf')]
        )
        self.onto = get_ontology(self.ontology_path).load()
        with self.onto:
            bandas_numericas = []
            for entidade in list(self.onto.individuals()) + list(self.onto.classes()):
                nome = entidade.name
                if nome.startswith("Banda"):
                    try:
                        valor = int(nome.replace("Banda_", "").replace("Banda", ""))
                        bandas_numericas.append(valor)
                    except ValueError:
                        continue
        # self.bandas_numericas = np.asarray(sorted(set(bandas_numericas)))
        self.bandas_numericas = np.array([ 400,  613,  622,  630, 1396, 1403, 1409, 1410, 1413, 1414, 1909,
                   1910, 1915, 1920, 1928, 1970, 2163, 2166, 2196, 2205, 2311,
                   2325, 2383, 2389])
    
    def registro(self, nome_amostra, nome_afloramento, bandas_detectadas, threshold=3):
        new_bands = []
        for band in bandas_detectadas:
            diff = abs(self.bandas_numericas - band)
            if diff.min() <= threshold:
                new_bands.append(self.bandas_numericas[np.argmin(diff)])
            else:
                new_bands.append(band)
        bandas_detectadas = new_bands
        
        onto = get_ontology(self.ontology_path).load()
        
        with onto:
            # 1. Localiza as classes principais
            classe_amostra = next((cls for cls in onto.classes() if cls.name == "Amostra"), None)
            classe_afloramento = next((cls for cls in onto.classes() if cls.name == "Afloramento"), None)
            Banda_Absorcao = next((cls for cls in onto.classes() if cls.name == "Banda_Absorcao"), None)
        
            # 2. Propriedades
            prop_tem_banda = onto.search_one(iri="*temBandaAbsorcao")
            prop_tem_banda.domain = [classe_amostra]
            prop_tem_banda.range = [Thing]
        
            prop_pertence_afloramento = onto.search_one(iri="*pertenceAoAfloramento")
            prop_pertence_afloramento.domain = [classe_amostra]
            prop_pertence_afloramento.range = [classe_afloramento]
        
            # 3. Cria ou localiza o indivíduo da amostra
            amostra_ind = onto.search_one(iri=f"*{nome_amostra}") or classe_amostra(nome_amostra)
        
            # 4. Cria ou localiza o afloramento
            afloramento_ind = onto.search_one(iri=f"*{nome_afloramento}") or classe_afloramento(nome_afloramento)
            amostra_ind.pertenceAoAfloramento = [afloramento_ind]
        
            # 5. Para cada banda
            for banda in bandas_detectadas:
                nome_classe_banda = f"Banda_{int(banda)}"
                nome_individuo_banda = f"Banda{int(banda)}"
        
                # 5.1 Cria a subclasse Banda_XXXX de Banda_Absorcao (caso não exista)
                classe_banda = onto.search_one(iri=f"*{nome_classe_banda}") or types.new_class(nome_classe_banda, (Banda_Absorcao,))
        
                # 5.2 Cria o indivíduo BandaXXXX como instância da subclasse Banda_XXXX
                banda_ind = onto.search_one(iri=f"*{nome_individuo_banda}") or classe_banda(nome_individuo_banda)
        
                # 5.3 Associa esse indivíduo à amostra
                if banda_ind not in getattr(amostra_ind, "temBandaAbsorcao", []):
                    amostra_ind.temBandaAbsorcao.append(banda_ind)


    def salvar_ontologia(self, caminho=None):
        print('Salvar_Ontologia')
        if not caminho:
            caminho = filedialog.asksaveasfilename(
                parent=None,
                title="Salvar Ontologia",
                defaultextension=".rdf",
                filetypes=[("Ontology files", "*.owl *.rdf")]
            )
        self.onto.save(file=caminho, format="rdfxml")








