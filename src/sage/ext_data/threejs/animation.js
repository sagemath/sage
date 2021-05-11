// animation.js

window.updateAnimation = (function() {

    var frameDuration = options.delay;
    var frameCount, totalDuration;
    var mixer, action, clock;
    var ui, slider;

    init();
    options.animationControls && initUI();
    options.autoPlay && play();
    return update;

    function init() {

        var frames = partitionScene();
        frameCount = frames.length;
        totalDuration = frameCount * frameDuration;

        var tracks = createTracks( frames );
        var clip = new THREE.AnimationClip( 'sage_animation', totalDuration, tracks );

        mixer = new THREE.AnimationMixer( scene );
        mixer.addEventListener( 'finished', function() { pause(); } );

        action = mixer.clipAction( clip );
        action.timeScale = 100; // Our times are in hundredths of a second.
        action.clampWhenFinished = true; // Sets paused=true instead of disabling.
        action.loop = options.loop ? THREE.LoopRepeat : THREE.LoopOnce;

        // We start/stop the clock to play/pause. The action's `paused` property
        // is instead used to determine if the animation has finished.
        action.play();

        clock = new THREE.Clock( false ); // Don't start playing yet.

    }

    // Group objects in by keyframe. The result is a potentially sparse list of
    // lists of objects. Assumes keyframe index specified in objects' user data.
    function partitionScene() {
        var frames = [];
        for ( var i=0 ; i < scene.children.length ; i++) {
            var obj = scene.children[i];
            var k = obj.userData && parseInt( obj.userData.keyframe );
            if ( k >= 0 ){
                if ( k in frames ) frames[k].push( obj );
                else frames[k] = [obj];
            }
        }
        return frames;
    }

    // Create and return a list of keyframe tracks to animate objects'
    // visibility depending on which keyframe they appear in.
    function createTracks( frames ) {
        var tracks = [];
        for ( var keyframe=0 ; keyframe < frames.length ; keyframe++ ) {
            var frame = frames[keyframe];
            if ( frame ) {
                // Show during the frame.
                var times = [keyframe * frameDuration];
                var values = [true];
                // For all but the first frame, start off hidden.
                if ( keyframe > 0 ) {
                    times.unshift( 0 );
                    values.unshift( false );
                }
                // For all but the last frame, end up hidden.
                if ( keyframe < frames.length - 1 ) {
                    times.push( ( keyframe + 1 ) * frameDuration );
                    values.push( false );
                }
                for ( var i=0; i < frame.length ; i++ ) {
                    var binding = frame[i].uuid + '.visible';
                    var track = new THREE.BooleanKeyframeTrack( binding, times, values );
                    tracks.push( track );
                }
            }
        }
        return tracks;
    }

    function initUI() {

        ui = document.getElementById( 'animation-ui' );

        initSlider();

        hookupButton( 'play-pause', togglePlayback );
        hookupButton( 'previous', gotoPreviousFrame );
        hookupButton( 'next', gotoNextFrame );
        hookupButton( 'slower', slowDown );
        hookupButton( 'faster', speedUp );
        hookupButton( 'toggle-loop', toggleLooping );

        function hookupButton( _class, onclick ) {
            var button = ui.getElementsByClassName( _class )[0];
            button.addEventListener( 'click', onclick );
            return button;
        }

    }

    function initSlider() {

        slider = ui.getElementsByClassName( 'slider' )[0];
        slider.value = action.time;
        slider.setAttribute( 'max', totalDuration );
        slider.addEventListener( 'change', change );
        slider.addEventListener( 'input', change );
        slider.addEventListener( 'mousedown', mousedown );
        slider.addEventListener( 'mouseup', mouseup );
        slider.addEventListener( 'keydown', keydown );

        function change() {
            // Keep animation in sync with position of the knob.
            seek( parseFloat( slider.value ) );
        }

        // Temporarily pause playback while the knob is being dragged around.
        var dragging = false;
        var wasPlaying = false;
        function mousedown( event ) {
            if ( event.button === 0 && !dragging ) {
                dragging = true;
                wasPlaying = isPlaying();
                pause();
            }
        }
        function mouseup( event ) {
            if ( event.button === 0 && dragging ) {
                dragging = false;
                if ( wasPlaying && !isFinished() ) play();
            }
        }

        function keydown( event ) {
            switch ( event.key ) {
                // Toggle playback using spacebar/enter.
                case ' ':
                case 'Enter':
                    if ( !event.repeat ) {
                        togglePlayback();
                        event.preventDefault();
                    }
                    break;
                // Adjust speed using up/down arrows.
                case 'Down':
                case 'ArrowDown':
                    slowDown();
                    event.preventDefault();
                    break;
                case 'Up':
                case 'ArrowUp':
                    speedUp();
                    event.preventDefault();
                    break;
                // Navigate frame-by-frame using left/right arrows.
                case 'Left':
                case 'ArrowLeft':
                    gotoPreviousFrame();
                    event.preventDefault();
                    break;
                case 'Right':
                case 'ArrowRight':
                    gotoNextFrame();
                    event.preventDefault();
                    break;
            }
        }

    }

    function isPlaying() {
        return clock.running;
    }

    function isFinished() {
        return action.paused;
    }

    function play() {
        if ( isFinished() ) action.reset(); // Return to the beginning.
        clock.start();
        updateUI();
        requestAnimationFrame( render ); // Re-enter render loop.
    }

    function pause() {
        clock.stop();
        updateUI();
    }

    function togglePlayback() {
        isPlaying() ? pause() : play();
    }

    function seek( time ) {

        // Ensure time is in range: [0, totalDuration].
        time = Math.max( 0, Math.min( totalDuration, time ) );
        action.time = time;

        // Make sure the animation is treated as finished if set to the very
        // end of the animation.
        action.paused = ( !isLooping() && time == totalDuration );

        updateUI();

        if ( !isPlaying() ) {
            // Ensure the animation gets re-rendered at the new position.
            requestAnimationFrame( render ); // Trigger a re-render.
        }

    }

    function moveKeyframes( delta ) {
        var keyframe = Math.floor( action.time / frameDuration );
        keyframe += delta;
        // Ensure keyframe is in range: [0, frameCount].
        keyframe = Math.max( 0, Math.min( frameCount, keyframe ) );
        seek( keyframe * frameDuration );
    }

    function gotoPreviousFrame() {
        pause();
        moveKeyframes( -1 );
    }

    function gotoNextFrame() {
        pause();
        moveKeyframes( +1 );
    }

    function speedUp() {
        action.timeScale *= 1.25;
    }

    function slowDown() {
        action.timeScale /= 1.25;
    }

    function isLooping() {
        return action.loop === THREE.LoopRepeat;
    }

    function setLooping( looping ) {
        action.loop = looping ? THREE.LoopRepeat : THREE.LoopOnce;
        updateUI();
    }

    function toggleLooping() {
        setLooping( !isLooping() );
    }

    function update() {
        mixer.update( clock.getDelta() );
        updateUI();
        return clock.running; // Should exit render loop if false.
    }

    function updateUI() {
        if ( ui ) {
            ui.classList.remove( 'playing', 'paused', 'once', 'loop' );
            ui.classList.add( isPlaying() ? 'playing' : 'paused' );
            ui.classList.add( isLooping() ? 'loop' : 'once' );
        }
        if ( slider ) slider.value = action.time;
    }

})();
