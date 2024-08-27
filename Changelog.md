# Changelog

## [1.0.0]alpha - 2024-08-26
No backward compatibility on nodal loads and boundary conditions definition.
- cambio externo: seteo de cargas y condiciones de vínculo: cambió el tipo de índice para la etiqueta local del GL - ya no hay compatibilidad hacia atrás en el 
- cambio externo: implementación de rutinas para hallar númedor de ecuaciones según nodo y GL local
- cambio externo: implementación de rutinas para hhallar desplazamientos y reacciones en GGL específicos
- cambio interno: cambio en la búsqueda de nodos por etiqueta externa:
  > implementación de findNode() en Struc2D3Dof
  > cambió create_elm
  > cambió create_NL
  > cambió create_BC

## [0.0.0]alpha - 2024-08-24
### Added
- First version with basic functionalities.
- Intended to model 2D framed structures with nodal loads.
- One simple example showing basic use.
