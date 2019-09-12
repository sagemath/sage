
==================================
Three.js JavaScript WebGL Renderer
==================================

A web-based interactive viewer using the Three.js JavaScript library maintained
by https://threejs.org.

The viewer is invoked by adding the keyword argument ``viewer='threejs'`` to the command 
``show()`` or any three-dimensional graphic. The scene is rendered and displayed
in the users's web browser. Interactivity includes

- Zooming in or out with the mouse wheel or pinching on a touch pad

- Rotation by clicking and dragging with the mouse or swiping on a touch pad

- Panning by right-clicking and dragging with the mouse or swiping with three fingers
  on a touch pad

The generated HTML file contains all data for the scene apart from the JavaScript library
and can be saved to disk for sharing or embedding in a web page. The option ``online``
can be set to ``true`` to provide links to the required files in an online content delivery
network. Alternately the required files can be downloaded from the Three.js GitHub repository
and linked directly from the web server.

Options currently supported by the viewer:

- ``aspect_ratio`` -- (default: [1,1,1]) list or tuple of three numeric
  values; `z`-aspect is automatically reduced when large but can be overridden

- ``axes`` -- (default: False) Boolean determining whether coordinate axes are drawn

- ``axes_labels`` -- (default: ['x','y','z']) list or tuple of three strings;
  set to False to remove all labels

- ``color`` -- (default: 'blue') color of the 3d object

- ``decimals`` -- (default: 2) integer determining decimals displayed in labels

- ``frame`` -- (default: True) Boolean determining whether frame is drawn

- ``online`` -- (default: False) Boolean determining whether the local standard package
  files are replaced by links to an online content delivery network

- ``opacity`` -- (default: 1) numeric value for transparency of lines and surfaces

- ``projection`` -- (default: 'perspective') the type of camera projection to use;
  'perspective' or 'orthographic'

- ``radius`` -- (default: None) numeric value for radius of lines; use to render
  lines thicker than available using ``thickness`` or on Windows platforms where
  ``thickness`` is ignored

- ``thickness`` -- (default: 1) numeric value for thickness of lines

AUTHORS:

- Paul Masson (2016): Initial version

EXAMPLES:

Three spheres of different color and opacity::

    sage: p1 = sphere(color='red', opacity='.5')
    sage: p2 = sphere((-1,-1,1), color='cyan', opacity='.3')
    sage: p3 = sphere((1,-1,-1), color='yellow', opacity='.7')
    sage: show(p1 + p2 + p3, viewer='threejs')

.. RAW:: html
    :file: threejs_examples/spheres.html

A parametric helix::

    sage: parametric_plot3d([cos(x),sin(x),x/10], (x,0,4*pi), color='red', viewer='threejs')
    Graphics3d Object

.. RAW:: html
    :file: threejs_examples/helix.html



.. RAW:: html

    <script>

    // iOS iframe auto-resize workaround

    if ( /(iPad|iPhone|iPod)/g.test( navigator.userAgent ) ) {

        var scenes = document.getElementsByTagName( 'iframe' );

        for ( var i=0 ; i < scenes.length ; i++ ) {

            scenes[i].style.width = getComputedStyle( scenes[i] ).width;
            scenes[i].style.height = getComputedStyle( scenes[i] ).height;
            scenes[i].setAttribute( 'scrolling', 'no' );

        }
    }

    </script>
