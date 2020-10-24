on run (volname, appname)

    tell application "Finder"
        tell disk (volname as string)
            open
                set opts to the icon view options of container window
                tell opts
                    set arrangement to not arranged
                    set icon size to 88
                end tell

                # Path .background/sage.png needs to be converted
                #  .background:sage.png
                set background picture of opts to file ".background:sage.png"

                tell container window
                    set current view to icon view
                    set toolbar visible to false
                    set statusbar visible to false
                    # The size of the window should match that of the image
                    set the bounds to {0, 0, 563, 348}

                    set position of item "README.txt" to {462, 52}
                    set position of item "Applications" to {462, 212}
                    set position of item (appname as string) to {104, 212}
                end tell

                # Some people report that the above settings are not applied
                # unless we wait for a couple of seconds in some versions of
                # Mac OS X.
                update without registering applications
                delay 3

            close
        end tell
    end tell

end run