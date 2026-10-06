# Análisis de Clústeres Biosintéticos de Tiopéptidos
**Autor:** Fernanda Rivas Romero  

---

## **METODOLOGÍA EMPLEADA HASTA AHORA**
### ***1. Búsqueda de tiopéptidos***
Dentro de la búsqueda de la literatura, la variedad de tiopéptidos es muy amplia dentro de los subgrupos denominados series. En este caso, la series d-e son las que se encuentran mayor caracterizadas, a diferencia de las series a, b y c. Es por ello que se realizó la búsqueda de tiopéptidos de todas las series para conocer que tan bien se encuentran caracterizados, esta caracterización enfocada en las seis enzimas indispensables para su macrociclización (YcaO, Ocin-Thiflike, FMN-dehydrogenase, dehydratase glutamylation and elimination domains, and [4 + 2]-macrocyclization), universalmente conservadas, a diferencia de las enzimas de modificación post-traduccional, las cuales varían significativamentre entre las series, e incluso dentro de los miembros de una misma serie. 

El criterio de selección de las seis enzimas indispensables para la macrociclización es que se encuentren registradas ya sea en MIBiG o en el NCBI con su número de GenBank (ACN...). En la búsqueda de artículos el nombre de la proteína codificada por un gen en específico para cada uno de los conjuntos de genes indispensables sí se encontraba, pero para otros sólo había el ID de la proteína (PYC...). Los tiopéptidos con los que se trabajarán para la búsqueda de homólogos están anotados en la tabla 1 (tiopeptidos-caracterizados.docx) con distintos metadatos: serie, organismo del cual fue aisalado, potencial antimicrobiano, otros roles además de su actividad antimicrobiana, ID en el repositorio MIBiG, y los genes biosintéticos para cada tiopéptido.

![Tabla 1. Tiopéptidos caracterizados que se emplearán para hacer los árboles filogenéticos](/Users/fernandarivasromero/Library/CloudStorage/OneDrive-Personal/Documentos/Maestría/UNAM/Documentos/avances/organismos/tabla1.png)
*Tabla 1. Metadatos de tiopéptidos de distintas series con sus genes biosintéticos indispensables para la macrociclización.*

Aunque las seis enzimas de la macrociclización son importantes para la formación de un tiopéptido, las enzimas [4+2] cicloadición y los dominios de eliminación y glutaminación de la deshidratasa son enzimas definitorias para la biosíntesis de tiopétidos.

De acuerdo con Just-Baringo et al., 2014 dentro de la biosíntesis del tiopéptido, la fosforilación y eliminación de Ser y Thr da lugar a los residuos de dehidroalanina (Dha) y dehidrobutirina (Dhb) respectivamente. La formación de Dha y Dhb permite que se lleve a cabo la cicloadición intramolecular del tipo aza-Diels-Alder entre residuos Dha, seguido de la deshidratación, y, cuando es necesario, de una eliminación para formar el anillo de seis miembros (Figura 1, Figura 2).

![Tiomuracina ejemplo 1](/Users/fernandarivasromero/Library/CloudStorage/OneDrive-Personal/Documentos/Maestría/UNAM/Documentos/avances/organismos/macrociclizacion.png)
*Figura 1. Ruta biosintética de la estructura básica de la tiomuracina que muestra el orden obligatorio de la formación de tiazoles, deshidratación y ciclización*

![Tiomuracina ejemplo 2](/Users/fernandarivasromero/Library/CloudStorage/OneDrive-Personal/Documentos/Maestría/UNAM/Documentos/avances/organismos/tiomuracina2.png)
*Figura 2. Otra ejemplificación de la ruta biosintética de la estructura básica de la tiomuracina que muestra la formación del anillo central a partir de la deshidratación y ciclización*

### ***2. BLAST***
De acuerdo con Aggarwal et al., 2021, la búsqueda de nuevos conjuntos de genes biosintéticos (CBGs) puede ser a través de la utilización de una secuencia concenso (query sequence) de la(s) proteína(s) que se considera(n) importante(s) para la biosíntesis de ese tiopéptido. En este caso, se tomaron las seis enzimas indispensables para su macrociclización. Aunque todas son importantes, en el artículo de Aggarwal et al., 2021 toman dos principales para la búsqueda de CGBs: **[4 + 2]-macrocyclization** y **dehydratase glutamylation and elimination domains**. La primera de ellas por ser quien le da la característica estructural del macrociclo, y es la enzima definitoria para la biosíntesis del tiopéptido. La segunda, por ser aquella que forma Dha (residuo de dehidroalanina), formada a partir de la deshidratación de residuos de serina, y que se requiere para la realizar la reacción formal aza-Diels-Alder para la formación del macrociclo.

La búsqueda de estas secuencias de proteínas homólogas en otros organismos se puede realizar a través de BLASTp (protein blast), donde los parámetros elegidos son:
1. **Choose Search Set:** ClusteredNR. Se elige esta base de datos que incluye un subset reducido pero representativo de la base de datos nr (NON-REDUNDANT). Ofrece una cobertura taxonómica más amplia ya que la búsqueda de secuencias lo hace también en homóloglos distantes. Además, su velocidad de búsqueda es más rápida que los non-redundant
2. **Program Selection:** PSI-BLAST. Se elige ya que realiza una búsqueda iterativa para una matriz de puntuación específica por posición (PSSM) para identificar parientes distantes en una familia de proteínas. Funciona si el tiopéptido es altamente divergente entre especies y blastp falla en la selección de homólogos distantes.

### ***3. Parámetros de filtrado y descarga de archivos fasta***
Para ahorrar tiempo de descarga, cuando PSI-BLAST arroja los resultados, en el apartado de "Filter Results" se hizo el filtrado para los parámetros de "Percent Identity" y "Query Coverage". Para el primero de 40-100, y para el segundo del 80-100. 

Una vez con el filtrado se pueden descargar todas las secuencias (select all) en formato FASTA (cluster) para obtener un formato multifasta de todas las secuencias a analizar. 

> *Nota: en BLAST hay Cluster Composition donde se engloban las secuencias de organismos homólogos a la proteína que estamos buscando, sin embargo, cuando se hace la descarga de los formatos FASTA, si el cluster esta compuesto de cuatro organismos, sólo se descarga uno, el más representativo, es decir, BLAST selecciona la secuencia más completa, con mejor cobertura o anotación dentro de este cluster. Esto evita que en el análisis MSA dominen secuencias casi idénticas (cepas o aislados sobre representados en la base de datos)*

### ***4. Generación de árboles filogenéticos de los tiopéptidos seleccionados***
Con los resultados del BLAST para las enzimas [4+2] cicloadición y los dominios de eliminación y glutaminación de la deshidratasa se realizaron árboles filogenéticos para ver la distribución los clados homólogos con las secuencias consenso. Se empleó el programa phylogeny.fr (https://phylogeny.fr/) en el apartado *"Phylogeny analysis: a la carte"* con los siguientes parámetros:

1. **Multiple Aligment:** se empleó MUSCLE debido a la precisión para conjuntos de datos pequeños-medianos (menos de 200-500 secuencias)
2. **Curación del alineamiento (opcional):** se empleó ClipKit debido a que es una herrmienta que preserva la precisión filogenética y la longitud de las ramas y no es tan agresivo como otro métodos de curación del alineamiento
3. **Árbol filogenético:** se empleó IQTree, un enfoque de Máxima Verosimilitud (ML) que busca la mejor combinación y selección automática de los modelos, velocidad y rigor estadístico. Además, implementa UFBoot (aproximación a bootstrap) que calcula el soporte de ramas de manera más rápida y con menor sesgo

### ***Anexo de conjuntos de genes biosintéticos (CGB)***
Este apartado es para poner las fotografías y referencias de donde se obtuvieron los nombres de las proteínas que conforman las seis enzimas indispensables para la macrociclización 

![CGBs de tiopeptidos de distintas series](/Users/fernandarivasromero/Library/CloudStorage/OneDrive-Personal/Documentos/Maestría/UNAM/Documentos/avances/organismos/thiopeptides-all.png)
*CGBs de tiopeptidos de distintas series con el nombre de sus proteínas que codifican para la macrociclización y post-traduccionales. Tomado de Vinogradov & Suga, 2020.*

![CGB tiostreptona.](/Users/fernandarivasromero/Library/CloudStorage/OneDrive-Personal/Documentos/Maestría/UNAM/Documentos/avances/organismos/Thriostepton-A.png)
*CGB tiostreptona. Tomado de Aggarwal et al., 2021.*

![CGB tiopeptina A1a.](/Users/fernandarivasromero/Library/CloudStorage/OneDrive-Personal/Documentos/Maestría/UNAM/Documentos/avances/organismos/thiopeptin-A1a.png)
*CGB tiopeptina A1a. (A) Distintas series de tiopéptidos clasificadas de acuerdo al cambio de oxidación en el anillo central. (B) Cambio de serie b a serie a con modificaciones en R y TpnL. (C) CGB de la triosteptona y tiopeptina indicando las posiciones de los genes homólogos con las líneas punteadas. Tomado de Ichikawa et al., 2018.*

![CGB kocurina y GE2270A.](/Users/fernandarivasromero/Library/CloudStorage/OneDrive-Personal/Documentos/Maestría/UNAM/Documentos/avances/organismos/kocurina_GE2270A.png)
*CGBs de la kocurina y GE2270A con sus genes homólogos. Tomado de Linares-Otoya et al., 2017.*
