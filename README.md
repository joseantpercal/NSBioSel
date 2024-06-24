# NSBioSel
Código fuente de experimentación para el TFM del Máster en Ciencia de Datos e Ingeniería de Computadores (UGR)

La implementación de los algoritmos utilizados en este estudio se ha basado en el repositorio disponible en GitHub, proporcionado por el Instituto J. Craig Venter. Los códigos fuente y las metodologías específicas pueden encontrarse en la siguiente dirección: https://github.com/JCVenterInstitute/NSForest/tree/master.

El presente trabajo se fundamenta en el artículo titulado "A machine learning method for the discovery of minimum marker gene combinations for cell type identification from single-cell RNA sequencing" de Aevermann et al. (2021). https://doi.org/10.1101/gr.275569.121.

Los datos empleados en este estudio han sido recopilados del artículo "A comprehensive mouse kidney atlas enables rare cell population characterization and robust marker discovery" de Novella-Rausell et al., 2023. https://doi.org/10.1016/j.isci.2023.106877.

## Autores

* José Antonio Pérez Calderón

## Tutores

* Luis Javier Herrera Maldonado
* Francisco M. Ortuño Muñoz

# Descripción

La identificación precisa de tipos celulares a partir de datos de secuenciación de ARN de célula única (scRNA-seq) es un desafío crucial en biología celular y medicina personalizada. El código contenido en este repositorio fue utilizado para realizar la experimentación del Trabajo fin de Máster. Este trabajo se centra en la evaluación y mejora del método NS-Forest, un algoritmo de aprendizaje automático conocido por su capacidad para seleccionar combinaciones mínimas de genes marcadores que optimizan la identificación de tipos celulares.

Para ello, se han implementado y evaluado algoritmos adicionales de selección de variables, como mRMR (minimum redundancy maximum relevance) y AUC-ROC, utilizando datos de test para evaluar su rendimiento en diferentes contextos. Esta evaluación incluyó analizar su efectividad en términos de precisión y sensibilidad para la discriminación de tipos celulares a través de diferentes conjuntos de datos. Se emplearon técnicas como "myDecisionTreeEvaluationTest" y "myDecisionTreeEvaluationTestMoved" para proporcionar una evaluación integral de la robustez y adaptabilidad de los modelos.
