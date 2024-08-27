# Changelog

## [1.0.0]alpha - 2024-08-27
No backward compatibility on nodal loads and boundary conditions definition.

Cambios "externos":
- seteo de cargas y condiciones de vínculo: cambió el tipo de índice para la etiqueta local del GL - ya no hay compatibilidad hacia atrás en el 
- implementación de métodos para retrieve desplazamientos y reacciones en GGL específicos: ask4U() y ask4R() en Struc2D3Dof
- implementación de método para retrieve elementos por etiqueta externa: Struc2D3Dof.ask4elem()
- implementación de método para retrieve carga nodal elemental por etiqueta externa de elemento, etiqueta local de nodo y etiqueta local de GL: Struc2D3Dof.ask4ENL()

Cambios "internos":
- cambio en la búsqueda de nodos por etiqueta externa:
  > implementación de findNode() en Struc2D3Dof
  > cambió create_elm
  > cambió create_NL
  > cambió create_BC
- implementación de método Struc2D3Dof.getEqLab() para hallar númedor de ecuaciones según nodo y GL local
- implementación de métodos getU() y getR() de Struc2D3Dof para hallar desplazamientos y reacciones según nodo y GL local
- implementación de método Struc2D3Dof.findElmts() para hallar elementos por etiqueta externa
- implementación de método ElmFrame2D.getL() para hallar carga nodal según nodo y GL local (en coordenadas globales y locales)


## [0.0.0]alpha - 2024-08-24
### Added
- First version with basic functionalities.
- Intended to model 2D framed structures with nodal loads.
- One simple example showing basic use.
