// fat_lines.js

var _fatLines = [];

function _createFatLine( lineStrip, geometry, materialOptions ) {

    var vertexCount = geometry.vertices.length;
    var positions = [];
    for ( var i=0 ; i < vertexCount ; i++ ) {
        var v = geometry.vertices[i];
        positions.push( v.x, v.y, v.z );
        // For a line strip, duplicate all but the first and last vertices.
        if ( lineStrip && i > 0 && i < vertexCount - 1 ) {
            positions.push( v.x, v.y, v.z );
        }
    }

    geometry = new THREE.LineSegmentsGeometry();
    geometry.setPositions( positions );

    var material = new THREE.LineMaterial( materialOptions );
    material.resolution = new THREE.Vector2( window.innerWidth, window.innerHeight );

    var line = new THREE.LineSegments2( geometry, material );
    line.computeLineDistances();
    line.scale.set( 1, 1, 1 );

    _fatLines.push( line );

    return line;

}

function createFatLineStrip( geometry, materialOptions ) {
    return _createFatLine( true, geometry, materialOptions );
}

function createFatLineSegments( geometry, materialOptions ) {
    return _createFatLine( false, geometry, materialOptions );
}

function rescaleFatLines() {
    var res = new THREE.Vector2( window.innerWidth, window.innerHeight );
    var n = _fatLines.length;
    for ( var i=0 ; i < n ; i++ ) {
        _fatLines[i].material.resolution = res;
    }
}
