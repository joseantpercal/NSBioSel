# NSBioSel
Código fuente de experimentación para el TFM del Máster en Ciencia de Datos e Ingeniería de Computadores (UGR)

## Autores

* José Antonio Pérez Calderón

## Tutores

* Luis Javier Herrera Maldonado
* Francisco M. Ortuño Muñoz

# Descripción

La identificación precisa de tipos celulares a partir de datos de secuenciación de ARN de célula única (scRNA-seq) es un desafío crucial en biología celular y medicina personalizada. El código contenido en este repositorio fue utilizado para realizar la experimentación del Trabajo fin de Máster. Este trabajo se centra en la evaluación y mejora del método NS-Forest, un algoritmo de aprendizaje automático conocido por su capacidad para seleccionar combinaciones mínimas de genes marcadores que optimizan la identificación de tipos celulares.

Para ello, se han implementado y evaluado algoritmos adicionales de selección de variables, como mRMR (minimum redundancy maximum relevance) y AUC-ROC, utilizando datos de test para evaluar su rendimiento en diferentes contextos. Esta evaluación incluyó analizar su efectividad en términos de precisión y sensibilidad para la discriminación de tipos celulares a través de diferentes conjuntos de datos. Se emplearon técnicas como "myDecisionTreeEvaluationTest" y "myDecisionTreeEvaluationTestMoved" para proporcionar una evaluación integral de la robustez y adaptabilidad de los modelos.
