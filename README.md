# Práctica: análisis filogenéticos con IQ-TREE y MrBayes

Esta práctica muestra un flujo de trabajo en terminal para:

1. organizar un proyecto de análisis filogenético;

2. analizar varios loci con IQ-TREE;

3. concatenar los loci con AMAS;

4. seleccionar un esquema/modelo con PartitionFinder;

5. configurar y ejecutar un análisis bayesiano con MrBayes.

Todos los pasos se ejecutan directamente en la terminal.
---

# 1. Descargar el epositorios

```bash
git clone https://github.com/cristoichkov/analisis_filogeneticos.git
cd analisis_filogeneticos
```

Para esta práctica agregaremos subdirectorios por etapa:

```bash
mkdir -p data/alacranes
mkdir -p out/01_iqtree
mkdir -p out/02_concatenacion
mkdir -p out/03_partitionfinder
mkdir -p out/04_mrbayes
```

La estructura será:

```text
analisis_filogeneticos/
├── bin/
├── contenedores 
├── data/
│   └── alacranes/
│       ├── 1410at6854.fasta
│       ├── 1941at6854.fasta
│       ├── 7960at6854.fasta
│       ├── 80357at6854.fasta
│       └── 9519at6854.fasta
├── out/
│   ├── 01_iqtree/
│   ├── 02_concatenacion/
│   ├── 03_partitionfinder/
│   └── 04_mrbayes/
└── README.md
```
---

# 2. Vincular los contenedores al proyecto

Para mantener el proyecto organizado sin duplicar archivos grandes, podemos crear un **enlace simbólico** hacia la carpeta donde ya están almacenados los contenedores de Apptainer.

Por ejemplo, si los contenedores se encuentran en:

```text
/srv/bishop/phylogenomics/contenedores/
```

desde la raíz de `analisis_filogeneticos` podemos crear:

```bash
ln -s /srv/bishop/phylogenomics/contenedores/ contenedores
```

Comprobar:

```bash
ls -lah contenedores/
```

Debería aparecer algo similar a:

```text
contenedores -> /srv/bishop/phylogenomics/contenedores/ 
```

Un **enlace simbólico** no copia los contenedores; solamente crea una referencia hacia la carpeta original.

Así podemos usar comandos cortos como:

```bash
apptainer exec contenedores/iqtree_3.1.3.sif iqtree3 --help
```

> Si el proyecto cambia de ubicación o se mueve a otra computadora, puede ser necesario crear nuevamente el enlace simbólico.

---

# 3. Revisar los alineamientos

Antes de hacer inferencia, conviene verificar cuántas secuencias contiene cada locus:

```bash
grep -c '^>' data/alacranes/*.fasta
```

Para revisar los nombres de las terminales:

```bash
grep '^>' data/alacranes/1410at6854.fasta
```

Todos los nombres dentro de un mismo FASTA deben ser únicos.

---

# 4. IQ-TREE: archivo de particiones

IQ-TREE puede analizar varios alineamientos como particiones independientes sin necesidad de concatenarlos previamente.

Crear:

```bash
nano meta/particiones_iqtree.nex
```

Contenido:

```nexus
#nexus
begin sets;
    charset locus_1410at6854  = data/alacranes/1410at6854.fasta: *;
    charset locus_1941at6854  = data/alacranes/1941at6854.fasta: *;
    charset locus_7960at6854  = data/alacranes/7960at6854.fasta: *;
    charset locus_80357at6854 = data/alacranes/80357at6854.fasta: *;
    charset locus_9519at6854  = data/alacranes/9519at6854.fasta: *;
end;
```

## ¿Qué significa `*`?

En:

```text
data/alacranes/1410at6854.fasta: *
```

el `*` significa:

> utilizar todas las posiciones del alineamiento.

Cada archivo FASTA se considera un locus completo.

---

# 5. Generar el archivo de particiones desde la terminal

Con cientos o miles de loci no conviene escribir cada línea manualmente.

Se puede generar el archivo directamente desde Bash:

```bash
echo '#nexus' > meta/particiones_iqtree.nex
echo 'begin sets;' >> meta/particiones_iqtree.nex
```

Después:

```bash
for f in data/alacranes/*.fasta
do
    locus=$(basename "$f" .fasta)
    echo "    charset locus_${locus} = ${f}: *;" >> meta/particiones_iqtree.nex
done
```

Finalmente:

```bash
echo 'end;' >> meta/particiones_iqtree.nex
```

Revisar:

```bash
cat meta/particiones_iqtree.nex
```

Contar las particiones:

```bash
grep -c 'charset' meta/particiones_iqtree.nex
```

---

# 6. Ejecutar IQ-TREE

Ejecutar:

```bash
apptainer exec contenedores/iqtree_3.1.3.sif \
    iqtree3 \
    -p meta/particiones_iqtree.nex \
    -m MFP \
    -B 1000 \
    --alrt 1000 \
    --partlh \
    -o Liocheles_australasiae \
    -T 4 \
    --prefix out/01_iqtree/alacranes
```

## Parámetros usados

| Parámetro | Significado |

|---|---|

| `-p` | Lee un archivo de particiones y utiliza el modelo particionado **edge-linked proportional**. |

| `-m MFP` | Ejecuta ModelFinder Plus para seleccionar automáticamente el modelo de sustitución. |

| `-B 1000` | Calcula 1000 réplicas de Ultrafast Bootstrap. |

| `--alrt 1000` | Calcula 1000 réplicas del SH-aLRT. |

| `--partlh` | Escribe la contribución de log-verosimilitud de cada partición. |

| `-o Liocheles_australasiae` | Define el outgroup con el que se escribe el árbol. |

| `-T 4` | Utiliza 4 hilos de CPU. |

| `--prefix` | Define la ruta y el prefijo de todos los archivos de salida. |

## UFBoot vs. bootstrap no paramétrico

En esta práctica utilizamos:

```bash
-B 1000
```

que corresponde a **Ultrafast Bootstrap (UFBoot)**. Lo elegimos porque es mucho más rápido y resulta adecuado para una sesión de clase.

IQ-TREE también permite realizar el **bootstrap no paramétrico clásico** con:

```bash
-b 1000
```

Por ejemplo:

```bash
apptainer exec contenedores/iqtree_3.1.3.sif \
    iqtree3 \
    -p meta/particiones_iqtree.nex \
    -m MFP \
    -b 1000 \
    --alrt 1000 \
    -o Liocheles_australasiae \
    -T 4 \
    --prefix out/01_iqtree/alacranes_bootstrap
```

La diferencia principal es:

```text
-B 1000  -> Ultrafast Bootstrap
-b 1000  -> bootstrap no paramétrico clásico
```

El bootstrap clásico suele tardar bastante más porque requiere búsquedas filogenéticas repetidas sobre matrices remuestreadas.

> **Para una publicación no es obligatorio utilizar bootstrap clásico.** UFBoot se utiliza ampliamente en trabajos publicados. Lo importante es reportar claramente el método empleado, el número de réplicas y, cuando sea posible, acompañarlo con otra medida de soporte como SH-aLRT. La elección entre UFBoot y bootstrap clásico depende del objetivo, tamaño del conjunto de datos y diseño del análisis.

En esta práctica usamos:

```text
SH-aLRT + UFBoot
```

---

## Salidas principales

```bash
ls out/01_iqtree/
```

Entre los archivos más importantes están:

```text
alacranes.treefile
alacranes.contree
alacranes.iqtree
alacranes.best_model.nex
alacranes.best_scheme
alacranes.partlh
alacranes.log
```

### Archivos de interés

- `alacranes.treefile`: árbol de máxima verosimilitud.

- `alacranes.contree`: árbol consenso con soportes.

- `alacranes.iqtree`: reporte detallado.

- `alacranes.best_model.nex`: modelos seleccionados.

- `alacranes.partlh`: log-verosimilitud por partición.

- `alacranes.log`: registro de la ejecución.

---

# 7. Concatenar los loci con AMAS

PartitionFinder necesita una matriz concatenada.

Para MrBayes también resulta conveniente trabajar con un único archivo NEXUS concatenado.

Generaremos:

```text
PHYLIP -> PartitionFinder
NEXUS  -> MrBayes
```

---

# 8. Crear el concatenado PHYLIP

Ejecutar:

```bash
apptainer exec contenedores/amas_1.0.sif \
    AMAS.py concat \
    -i data/alacranes/*.fasta \
    -f fasta \
    -d dna \
    -t out/02_concatenacion/alacranes_concat.phy \
    -p out/02_concatenacion/particiones_concat.txt \
    -u phylip
```

## Parámetros de AMAS

| Parámetro | Significado |

|---|---|

| `concat` | Concatena varios alineamientos. |

| `-i` | Especifica los archivos de entrada. |

| `-f fasta` | Indica que los alineamientos de entrada están en FASTA. |

| `-d dna` | Define el tipo de datos como DNA. |

| `-t` | Define el archivo concatenado de salida. |

| `-p` | Define el archivo donde AMAS escribirá las coordenadas de las particiones. |

| `-u phylip` | Escribe la matriz concatenada en formato PHYLIP. |

Revisar:

```bash
cat out/02_concatenacion/particiones_concat.txt
```

Para estos datos:

```text
p1_1410at6854 = 1-921
p2_1941at6854 = 922-1986
p3_7960at6854 = 1987-2730
p4_80357at6854 = 2731-3669
p5_9519at6854 = 3670-4446
```

La longitud total del superalineamiento es:

```text
4446 posiciones
```

---

# 9. Crear el concatenado NEXUS para MrBayes

Ejecutar:

```bash
apptainer exec contenedores/amas_1.0.sif \
    AMAS.py concat \
    -i data/alacranes/*.fasta \
    -f fasta \
    -d dna \
    -t out/02_concatenacion/alacranes_concat.nex \
    -p out/02_concatenacion/particiones_concat.txt \
    -u nexus
```

Ahora:

```bash
ls out/02_concatenacion/
```

debería contener:

```text
alacranes_concat.nex
alacranes_concat.phy
particiones_concat.txt
```

---

# 10. Preparar PartitionFinder

Crear una copia del alineamiento PHYLIP dentro de la etapa de PartitionFinder:

```bash
cp out/02_concatenacion/alacranes_concat.phy \
   out/03_partitionfinder/
```

Crear:

```bash
nano out/03_partitionfinder/partition_finder.cfg
```

Contenido:

```ini
alignment = alacranes_concat.phy;
branchlengths = linked;
models = mrbayes;
model_selection = BIC;
[data_blocks]
p1_1410at6854 = 1-921;
p2_1941at6854 = 922-1986;
p3_7960at6854 = 1987-2730;
p4_80357at6854 = 2731-3669;
p5_9519at6854 = 3670-4446;
[schemes]
search = greedy;
```

---

# 11. Parámetros de PartitionFinder

## `alignment`

```ini
alignment = alacranes_concat.phy;
```

Especifica la matriz concatenada.

## `branchlengths = linked`

```ini
branchlengths = linked;
```

Las particiones comparten las longitudes relativas de las ramas.

## `models = mrbayes`

```ini
models = mrbayes;
```

Restringe la selección a modelos que pueden implementarse posteriormente en MrBayes.

## `model_selection = BIC`

```ini
model_selection = BIC;
```

Utiliza el Bayesian Information Criterion.

El BIC considera:

- qué tan bien explica el modelo los datos;

- el número de parámetros del modelo;

- el tamaño de la matriz.

Por tanto, no siempre selecciona el modelo más complejo.

## `[data_blocks]`

Define los límites de cada locus dentro de la matriz concatenada.

## `search = greedy`

```ini
search = greedy;
```

Permite que PartitionFinder combine diferentes bloques si compartir un modelo produce un mejor BIC.

Por eso cinco loci iniciales pueden terminar agrupados en menos particiones.

---

# 12. Ejecutar PartitionFinder

Si existe una corrida anterior incompleta:

```bash
rm -rf out/03_partitionfinder/analysis
```

Ejecutar con cuatro procesos:

```bash
HDF5_USE_FILE_LOCKING=FALSE \
apptainer exec contenedores/partitionfinder2_2.1.1.sif \
    PartitionFinder.py \
    --processes 4 \
    out/03_partitionfinder
```

## Parámetros

| Parámetro | Significado |

|---|---|

| `HDF5_USE_FILE_LOCKING=FALSE` | Evita problemas de bloqueo de la base HDF5 en algunos sistemas de archivos. |

| `--processes 4` | Limita PartitionFinder a cuatro procesos de CPU. |

| `out/03_partitionfinder` | Directorio que contiene el `.cfg` y el alineamiento. |

---

# 13. Salida de PartitionFinder para MrBayes

En este ejemplo, PartitionFinder agrupó los cinco loci en un único subconjunto:

```nexus
begin mrbayes;
    charset Subset1 = 1-921 922-1986 1987-2730 2731-3669 3670-4446;
    partition PartitionFinder = 1:Subset1;
    set partition=PartitionFinder;
    lset applyto=(1) nst=6 rates=invgamma;
end;
```

Esto indica que PartitionFinder decidió que todos los loci podían analizarse con el mismo modelo.

La instrucción:

```nexus
nst=6
```

corresponde a un modelo GTR.

La instrucción:

```nexus
rates=invgamma
```

incorpora:

```text
+I + G
```

Por tanto, el modelo corresponde a:

```text
GTR + I + G
```

---

# 14. Preparar el archivo de MrBayes

Copiar el NEXUS concatenado a la etapa de MrBayes:

```bash
cp out/02_concatenacion/alacranes_concat.nex \
   out/04_mrbayes/
```

Editar:

```bash
nano out/04_mrbayes/alacranes_concat.nex
```

Al final del archivo, después del bloque `DATA`, agregar:

```nexus
begin mrbayes;
    set autoclose=yes nowarn=yes;
    charset Subset1 = 1-921 922-1986 1987-2730 2731-3669 3670-4446;
    partition PartitionFinder = 1:Subset1;
    set partition=PartitionFinder;
    lset applyto=(1) nst=6 rates=invgamma;
    outgroup Liocheles_australasiae;
    mcmcp ngen=1000000
          nruns=2
          nchains=4
          samplefreq=100
          printfreq=1000
          diagnfreq=5000
          burninfrac=0.25;
    mcmc;
    sump;
    sumt;
end;
```

---

# 15. Parámetros de MrBayes

## `set autoclose=yes`

```nexus
set autoclose=yes nowarn=yes;
```

Permite que MrBayes termine automáticamente después de ejecutar el archivo.

`nowarn=yes` reduce preguntas interactivas durante la ejecución.

---

## `charset`

```nexus
charset Subset1 = 1-921 922-1986 1987-2730 2731-3669 3670-4446;
```

Define qué posiciones pertenecen al subconjunto seleccionado por PartitionFinder.

---

## `partition`

```nexus
partition PartitionFinder = 1:Subset1;
set partition=PartitionFinder;
```

Crea el esquema de particionamiento y lo activa.

---

## Modelo de sustitución

```nexus
lset applyto=(1) nst=6 rates=invgamma;
```

Significa:

```text
nst=6          -> GTR
rates=invgamma -> sitios invariantes + distribución gamma
```

---

## Outgroup

```nexus
outgroup Liocheles_australasiae;
```

Define el grupo externo utilizado para presentar el árbol enraizado.

---

# 16. Configuración de la MCMC

```nexus
mcmcp ngen=1000000
      nruns=2
      nchains=4
      samplefreq=100
      printfreq=1000
      diagnfreq=5000
      burninfrac=0.25;
```

## `ngen=1000000`

Ejecuta:

```text
1,000,000 generaciones
```

---

## `nruns=2`

Ejecuta ****dos corridas MCMC independientes****.

Las dos corridas:

- parten de estados independientes;

- no intercambian estados entre sí;

- permiten comparar si ambas llegan a la misma distribución posterior.

---

## `nchains=4`

Cada corrida contiene cuatro cadenas:

```text
Corrida 1
├── 1 cadena fría
├── 3 cadenas calientes
│
Corrida 2
├── 1 cadena fría
└── 3 cadenas calientes
```

Por tanto:

```text
2 corridas × 4 cadenas = 8 cadenas
```

Las cuatro cadenas dentro de una misma corrida ****no son independientes****.

Las cadenas calientes y la cadena fría pueden intentar intercambiar estados.

Las cadenas calientes facilitan la exploración del espacio de:

- topologías;

- longitudes de rama;

- parámetros del modelo.

La cadena fría es la que muestrea la distribución posterior objetivo.

Las dos corridas (`nruns=2`), en cambio, sí son independientes una de otra.

---

## `samplefreq=100`

```nexus
samplefreq=100
```

Guarda una muestra cada 100 generaciones.

Con un millón de generaciones:

```text
1,000,000 / 100 = 10,000 muestras por corrida
```

antes del burn-in.

---

## `printfreq=1000`

Muestra el progreso en la terminal cada 1000 generaciones.

No cambia el análisis.

---

## `diagnfreq=5000`

Calcula diagnósticos de convergencia cada 5000 generaciones.

Con dos corridas, MrBayes puede calcular medidas como el:

```text
Average Standard Deviation of Split Frequencies
```

---

## `burninfrac=0.25`

Descarta el primer 25 % de las muestras antes de resumir la distribución posterior.

La idea es eliminar la fase inicial previa a la estabilización de la cadena.

---

# 17. Ejecutar MrBayes

Entrar al directorio de resultados:

```bash
cd out/04_mrbayes
```

Ejecutar con MPI.

Como existen:

```text
2 corridas × 4 cadenas = 8 cadenas
```

```bash
mpirun -np 4 \
    apptainer exec ../../contenedores/mrbayes_*.sif \
    mb alacranes_concat.nex
```

> `nchains=4` y `-np 4` no significan lo mismo.

>

> `nchains` controla el número de cadenas MCMC por corrida.

>

> `-np` controla el número de procesos MPI disponibles para ejecutar el análisis.

---

# 18. Archivos de salida de MrBayes

Después del análisis pueden aparecer archivos como:

```text
alacranes_concat.nex.p
alacranes_concat.nex.t
alacranes_concat.nex.con.tre
alacranes_concat.nex.parts
alacranes_concat.nex.trprobs
alacranes_concat.nex.pstat
alacranes_concat.nex.tstat
alacranes_concat.nex.vstat
alacranes_concat.nex.lstat
alacranes_concat.nex.mcmc
```

Los más importantes son:

- `.p`: muestras de los parámetros de la MCMC.

- `.t`: árboles muestreados.

- `.con.tre`: árbol consenso.

- `.trprobs`: probabilidades posteriores de árboles.

- `.pstat`: resumen de parámetros.

- `.tstat`: resumen de árboles.

---

# 19. Flujo completo

```text
data/alacranes/*.fasta
        │
        ├──────────────────────────────┐
        │                              │
        ▼                              ▼
meta/particiones_iqtree.nex           AMAS
        │                              │
        ▼                     ┌────────┴────────┐
     IQ-TREE                  ▼                 ▼
        │             alacranes_concat.phy   alacranes_concat.nex
        │                     │                 │
        ▼                     ▼                 │
   árbol ML             PartitionFinder         │
   ModelFinder                │                  │
   SH-aLRT                    ▼                  │
   UFBoot             modelo / esquema          │
   logL por locus             │                  │
                              └────────┬─────────┘
                                       ▼
                                    MrBayes
                                       │
                                       ▼
                          distribución posterior
                          probabilidades posteriores
                          árbol consenso
```

---

# 20. Estructura final del proyecto

```text
analisis_filogeneticos/
├── bin/
│   └── reporte.Rmd
├── contenedores -> /ruta/a/contenedores/
├── data/
│   └── alacranes/
│       ├── 1410at6854.fasta
│       ├── 1941at6854.fasta
│       ├── 7960at6854.fasta
│       ├── 80357at6854.fasta
│       └── 9519at6854.fasta
│
├── meta/
│   └── particiones_iqtree.nex
│
├── out/
│   ├── 01_iqtree/
│   │   ├── alacranes.treefile
│   │   ├── alacranes.contree
│   │   ├── alacranes.iqtree
│   │   └── ...
│   │
│   ├── 02_concatenacion/
│   │   ├── alacranes_concat.phy
│   │   ├── alacranes_concat.nex
│   │   └── particiones_concat.txt
│   │
│   ├── 03_partitionfinder/
│   │   ├── alacranes_concat.phy
│   │   ├── partition_finder.cfg
│   │   └── analysis/
│   │
│   └── 04_mrbayes/
│       ├── alacranes_concat.nex
│       ├── alacranes_concat.nex.p
│       ├── alacranes_concat.nex.t
│       ├── alacranes_concat.nex.con.tre
│       └── ...
│
└── README.md
```

---

# 21. Visualización y manipulación del árbol de IQ-TREE en R

La visualización en R se realizará utilizando el árbol obtenido con **IQ-TREE**, específicamente:

```text
out/01_iqtree/alacranes.contree
```

Para esta práctica, las instrucciones detalladas de lectura, manipulación y visualización del árbol con **R**, `treeio` y `ggtree` se encuentran en el archivo:

```text
bin/reporte.Rmd
```

Este R Markdown contiene los pasos para:

- leer el árbol de IQ-TREE con `read.iqtree()`;
- explorar la estructura del objeto `treedata`;
- visualizar la topología y las longitudes de rama;
- mostrar los soportes `SH-aLRT` y `UFBoot`;
- representar los soportes mediante símbolos;
- asociar información taxonómica a las terminales;
- colorear las terminales por familia;
- ajustar el lienzo y preparar figuras.

Para abrirlo:

```bash
ls bin/
```

y posteriormente trabajar con:

```text
bin/reporte.Rmd
```

---

# 22. Flujo final de la práctica

```text
Alineamientos FASTA
        │
        ├── IQ-TREE
        │     ├── ModelFinder
        │     ├── Máxima Verosimilitud
        │     ├── SH-aLRT
        │     └── UFBoot
        │
        └── AMAS
              │
              ├── PHYLIP
              │     └── PartitionFinder
              │             └── esquema/modelo
              │
              └── NEXUS
                    └── MrBayes
                          ├── MCMC
                          ├── distribución posterior
                          └── árbol consenso
Visualización:
out/01_iqtree/alacranes.contree
        │
        ▼
bin/reporte.Rmd
        │
        ▼
R + treeio + ggtree
```
