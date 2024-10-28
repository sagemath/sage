.. _threejs_viewer:

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
  values

- ``auto_scaling`` -- (default: [False, False, False]) list or tuple of three booleans;
  set to True to automatically scale down the corresponding direction if it is too
  large

- ``axes`` -- (default: ``False``) Boolean determining whether coordinate axes are drawn

- ``axes_labels`` -- (default: ['x','y','z']) list or tuple of three strings;
  set to False to remove all labels

- ``axes_labels_style`` -- (default: None) list of three dicts, one per axis, or
  a single dict controlling all three axes; supports the same styling options as
  :func:`~sage.plot.plot3d.shapes2.text3d` such as ``color``, ``opacity``, ``fontsize``,
  ``fontweight``, ``fontstyle``, and ``fontfamily``

- ``color`` -- (default: 'blue') color of the 3D object

- ``decimals`` -- (default: 2) integer determining decimals displayed in labels

- ``depth_write`` -- (default: ``True`` for opaque surfaces, False for transparent surfaces)
  whether to write the surface's depth into the depth buffer for the purpose of occluding
  objects behind it

- ``frame`` -- (default: ``True``) Boolean determining whether frame is drawn

- ``online`` -- (default: ``False``) Boolean determining whether the local standard package
  files are replaced by links to an online content delivery network

- ``opacity`` -- (default: 1) numeric value for transparency of lines and surfaces

- ``page_title`` -- (default: None) string containing the title of the generated HTML page; often
  displayed in the browser window's title bar, its tab list, and/or the operating system's task bar

- ``projection`` -- (default: 'perspective') the type of camera projection to use;
  'perspective' or 'orthographic'

- ``radius`` -- (default: None) numeric value for radius of lines; use to render
  lines thicker than available using ``thickness`` or on Windows platforms where
  ``thickness`` is ignored

- ``render_order`` -- (default: 0) numeric value for rendering order of transparent surfaces;
  objects render from lowest to highest value ensuring that lower-valued objects render completely

- ``single_side`` -- (default: ``False``) Boolean determining whether both sides of a surface material
  are rendered; set to True to reduce rendering artifacts for closed transparent surfaces

- ``theme`` -- (default: 'light') the color scheme to use for the scene and user interface;
  'light' or 'dark'

- ``thickness`` -- (default: 1) numeric value for thickness of lines

- ``viewpoint`` -- (default: None) list or tuple of the form [[x,y,z],angle] setting the initial
  viewpoint of the scene, where angle is in degrees; can be determined using the 'Get Viewpoint'
  option of the information menu

In addition, the following animation-related options are supported:

- ``animate`` -- (default: depends) whether to enable animation. Automatically set to ``True``
  if animation data is present in the plot. If ``False``, all frames of animation will be displayed
  simultaneously.

- ``animation_controls`` -- (default: ``True``) whether to include the playback slider and buttons
  (play, pause, etc.) in the page

- ``auto_play`` -- (default: ``True``) whether to immediately start playing the animation when the page
  loads. Recommend setting ``animation_controls=True`` to be able to start playback.

- ``delay`` -- (default: 20) an integer amount of time between consecutive frames of animation,
  in hundredths of a second

- ``loop`` -- (default: ``True``) whether to loop the animation or have it stop after reaching the end.
  Can be toggled on the page itself if ``animation_controls`` is set.

Clicking on the information icon in the lower right-hand corner of the viewer opens
a menu of available actions. These include saving the three-dimensional scene as a static
PNG image or as complete HTML source code.

AUTHORS:

- Paul Masson (2016): Initial version
- Joshua Campbell (2020): Animation support

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

An :meth:`~sage.plot.animate.Animation.interactive` animation::

  sage: def build_frame(t):
  ....:     e = parametric_plot3d([sin(x-t), 0, x], (x, 0, 2*pi), color='red')
  ....:     m = parametric_plot3d([0, -sin(x-t), x], (x, 0, 2*pi), color='green')
  ....:     return e + m
  sage: frames = [build_frame(t) for t in (0, pi/32, pi/16, .., 2*pi)]
  sage: plot = animate(frames).interactive()
  sage: show(plot, delay=5, auto_play=False, projection='orthographic')

.. RAW:: html
    :file: threejs_examples/animation.html


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
