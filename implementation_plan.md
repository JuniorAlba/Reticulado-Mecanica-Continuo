# Plan de Finalización — TP1 Mecánica del Continuo (Caso 12)

## Contexto

El proyecto simula la respuesta dinámica de un **reticulado plano 2D** (Caso 12, 15 nodos, 26 barras) con condiciones de empotramiento en nodos 1 y 15, bajo carga uniforme P = 3.2 aplicada en nodos 7 y 10.

### Archivos existentes

| Archivo | Estado | Descripción |
|---|---|---|
| `TP/caso12_inicio.m` | ✅ Completo | Geometría, rigidez, masas, frecuencias naturales |
| `TP/caso12_dinamica.m` | ⚠️ Incompleto / con bugs | Simulación dinámica lineal y no lineal |
| `Consigna/gif.m` | ✅ Provisto | Función para generar GIF animado |
| `Consigna/EjemploGif.m` | ✅ Provisto | Ejemplo de uso de gif.m |
| `Consigna/f_LectDxf.m` | ✅ Provisto | Lector de archivos DXF |

---

## Problemas Identificados en el Código Actual

### 1. Geometría inconsistente entre archivos
- `caso12_inicio.m` usa coordenadas `coords` (comentadas con los 15 nodos correctos)
- `caso12_dinamica.m` usa coordenadas `x0/y0` **diferentes** — los valores de y0 no coinciden (e.g. nodo 3 → y=5 en dinámica vs y=15 en inicio)

> [!CAUTION]
> La geometría en `caso12_dinamica.m` (x0/y0 en líneas 30–31) debe verificarse contra la consigna del PDF. Las coordenadas en `caso12_inicio.m` parecen más cuidadosas (con comentarios de nodos). Hay que unificarlas.

### 2. Detección de `tF` — Criterio de colapso
La función `detectar_tF` usa cambio de signo en el área de triángulos para detectar el colapso topológico. La lista de triángulos definida (líneas 59–74) podría estar incompleta o mal orientada para la geometría real del Caso 12.

### 3. Amortiguador (punto c) — Area de referencia
La fórmula del área de referencia `AR` en `drag_force` usa `pi*(r^2 - dist^2)` para la transición parcial. Esto no es la fórmula correcta para la sección transversal de una esfera cortada a una altura:
- La fórmula correcta para la sección circular a profundidad `d` desde el centro es: `AR = pi * (r^2 - d^2)` cuando `|d| < r` → ✅ esto **sí es correcto** si `dist` es la distancia del centro a la superficie.

> [!NOTE]
> Hay que revisar la convención de signos: `dist = h_SL - y_c`. Si `y_c < h_SL` → esfera bajo la superficie libre → sumergida. La lógica actual parece correcta pero requiere validación.

### 4. Animación — Falta integración con `gif.m`
La animación actual (líneas 307–332) usa `drawnow + pause` en lugar de guardar un GIF usando `gif.m` del profesor. El enunciado muy probablemente requiere entregar un GIF.

### 5. Faltan gráficas de frecuencias naturales
`caso12_inicio.m` calcula frecuencias pero no las grafica. El TP suele pedir el gráfico de modos propios.

### 6. No hay integración con DXF
La consigna provee `f_LectDxf.m` para leer la geometría desde un archivo `.dxf`. Hay que verificar si el TP requiere leer el reticulado desde el DXF en lugar de hardcodear las coordenadas.

---

## Open Questions

> [!IMPORTANT]
> **Geometría desde DXF**: ¿El TP requiere leer la geometría (nodos y barras) desde el archivo `.dxf` del Caso 12 usando `f_LectDxf.m`? ¿Tenés el archivo DXF del Caso 12?

> [!IMPORTANT]
> **Coordenadas correctas**: ¿Las coordenadas en `caso12_inicio.m` (nodo 3 en y=15, nodo 5 en y=15, nodo 8 en y=30) son las correctas según el PDF de la consigna? O ¿son las de `caso12_dinamica.m` (nodo 3 en y=5, etc.)?

> [!IMPORTANT]
> **Puntos del TP pendientes**: ¿Cuáles de estos puntos del enunciado ya están resueltos y cuáles faltan?
> - a.i) Hipótesis no lineal (grandes desplazamientos)
> - a.ii) Hipótesis lineal (pequeños desplazamientos)
> - b) Tensión en barra `a` y coordenada actual de nodo `b`
> - c) Amortiguador por fuerza de arrastre (Drag Force)
> - d) Animación GIF de la evolución temporal

---

## Tareas Propuestas

### FASE 1 — Corrección y unificación de la geometría

#### [MODIFY] [caso12_inicio.m](file:///c:/Users/joack/Documents/Reticulado-Mecanica-Continuo/TP/caso12_inicio.m)
- Verificar coordenadas contra el PDF
- Agregar gráfica de la estructura inicial con numeración de nodos y barras

#### [MODIFY] [caso12_dinamica.m](file:///c:/Users/joack/Documents/Reticulado-Mecanica-Continuo/TP/caso12_dinamica.m)
- Unificar coordenadas x0/y0 con las de `caso12_inicio.m`
- Verificar y corregir la lista de triángulos para detección de colapso

---

### FASE 2 — Completar la dinámica (puntos a.i, a.ii, b)

#### [MODIFY] [caso12_dinamica.m](file:///c:/Users/joack/Documents/Reticulado-Mecanica-Continuo/TP/caso12_dinamica.m)
- Verificar que `ode_lineal` y `ode_nolineal` estén correctamente formuladas
- Ajustar la detección de `tF` con triángulos válidos para Caso 12
- Agregar tabla comparativa: `tF_lin` vs `tF_nl`
- Agregar gráfica: tensión normal en barra `a` (vs tiempo) para las 3 simulaciones
- Agregar gráfica: coordenada actual `y` del nodo `b` (vs tiempo)

---

### FASE 3 — Amortiguador (punto c)

#### [MODIFY] [caso12_dinamica.m](file:///c:/Users/joack/Documents/Reticulado-Mecanica-Continuo/TP/caso12_dinamica.m)
- Revisar y documentar la función `drag_force`
- Verificar que el nodo del amortiguador (`nodo_c = 8`) sea el correcto
- Validar el efecto del amortiguador comparando `tF_am` vs `tF_nl`

---

### FASE 4 — Animación GIF (punto d)

#### [MODIFY] [caso12_dinamica.m](file:///c:/Users/joack/Documents/Reticulado-Mecanica-Continuo/TP/caso12_dinamica.m)
- Reemplazar la animación con `pause` por generación de GIF usando `gif.m`
- Generar GIF para la simulación no lineal (a.i)
- Copiar `gif.m` a la carpeta `TP/` o ajustar el path

#### [NEW] `TP/gif.m` (copia/symlink)
- Copiar `Consigna/gif.m` a `TP/` para que esté accesible desde el script principal

---

### FASE 5 — Lectura desde DXF (opcional / si la consigna lo requiere)

#### [NEW] `TP/caso12_desde_dxf.m`
- Leer la geometría del Caso 12 desde el archivo `.dxf` usando `f_LectDxf.m`
- Extraer nodos y conectividad automáticamente

---

### FASE 6 — Organización y entrega

#### [NEW] `TP/caso12_main.m` (script maestro)
- Un único script que llame a inicio + dinámica en orden
- Genere todos los resultados requeridos
- Guarde figuras como `.png` y GIF como `.gif`

#### [MODIFY] `README.md`
- Documentar la estructura del proyecto
- Instrucciones de ejecución

---

## Plan de Verificación

### Verificación automática (MATLAB)
1. Ejecutar `caso12_inicio.m` → verificar frecuencias naturales (6 modos)
2. Ejecutar `caso12_dinamica.m` → verificar que las 3 simulaciones corran sin error
3. Verificar que `tF_nl < tF_lin` (lo no lineal colapsa antes que lo lineal, generalmente)
4. Verificar que el GIF se genere correctamente

### Verificación manual
- Comparar gráficas de tensión y desplazamiento con resultados esperados
- Verificar visualmente la animación del reticulado deformándose
- Confirmar que el amortiguador reduce la amplitud o retrasa el colapso
