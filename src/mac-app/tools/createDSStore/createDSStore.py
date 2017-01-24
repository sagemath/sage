#!/usr/bin/python

"""
Creates a .DS_Store file that instructs Mac OS X's finder to use icon view with
a certain background image, and places the icons.

This file is intended to be placed into the directory from which a dmg is
created for packaging an application. It needs as additional information the
name of the volume that will be used when creating the dmg and the name of the
application, including the .app, e.g., Sage-6.9.app.

We could just always use the same .DS_Store that if the volume name and
application name were not changing.

Most of the other values set in the .DS_Store have been taken from dmgbuild
(on pypi). However, dmgbuild requires the creation of a writable dmg, then
puts the .DS_Store file into it, and converts it to a read-only dmg. We avoid
this by using the volume name of the dmg instead.
See getBackgroundImage_alias.
"""
from __future__ import print_function

import ds_store, mac_alias, biplist
import datetime, sys, os

# Hard-coded data for finder
# The .background.png's dpi should be set to 72, otherwise, we get weird
# behavior (especially on Retina displays). This can be fixed with
#    convert -density 72 -units pixelsperinch file.png file.png
backgroundImageSize = (563, 348)
backgroundImageName = '.background.png'
iconPositions = {
    'README.txt' : (497, 52),
    'Applications' : (462, 212),
    'app' :  (104, 212)
}

def getWindowData_bwsp():
    global backgroundImageSize

    return {
        'ShowTabView': False,
        'ShowStatusBar': False,
        'WindowBounds': '{{100, 100}, {%d, %d}}' % backgroundImageSize,
        'PreviewPaneVisibility': False,
        'ContainerShowSidebar': False,
        'ShowToolbar': False,
        'ShowPathbar': False,
        'SidebarWidth': 180,
        'ShowSidebar': False
    }

def getSomeTime():
    class UTC(datetime.tzinfo):
        def utcoffset(self, dt):
            return datetime.timedelta(0)
        def tzname(self, dt):
            return "UTC"
        def dst(self, dt):
            return datetime.timedelta(0)

    return datetime.datetime(2000, 1, 1, tzinfo = UTC())

def getBackgroundImage_alias(volume_name):
    """
    The location of the background image is stored as "Alias" which consists
    of the volume name and an absolute path with respect to that volume.
    
    Currently, we only support having the background image at the root of the
    volume.
    """

    # Also see Alias.for_file which only works after the volume is mounted

    global backgroundImageName

    carbon_path = '%s:%s' % (volume_name, backgroundImageName)

    volume = mac_alias.VolumeInfo(volume_name,
                                  getSomeTime(),
                                  'H+',
                                  0,
                                  0,
                                  '\x00\x00')

    target = mac_alias.TargetInfo(0,
                                  backgroundImageName,
                                  0, 0,
                                  getSomeTime(),
                                  '\x00\x00\x00\x00',
                                  '\x00\x00\x00\x00',
                                  folder_name = volume_name,
                                  carbon_path = carbon_path)

    return mac_alias.Alias(volume = volume, target = target)

def getBackgroundImage_icvp(volume_name):
    return {
        'gridSpacing': 100.0,
        'textSize': 16.0,
        'viewOptionsVersion': 1,
        'backgroundColorBlue': 1.0,
        'scrollPositionX': 0.0,
        'iconSize': 88.0,
        'backgroundColorGreen': 1.0,
        'arrangeBy': 'none',
        'showIconPreview': False,
        'gridOffsetX': 0.0,
        'gridOffsetY': 0.0,
        'showItemInfo': False,
        'labelOnBottom': True,
        'backgroundType': 2,
        'scrollPositionY': 0.0,
        'backgroundColorRed': 1.0,
        'backgroundImageAlias': biplist.Data(
            getBackgroundImage_alias(volume_name).to_bytes())
        }

def createDSStore(target_dir, volume_name, app_name):
    global iconPositions

    filePath = os.path.join(target_dir, '.DS_Store')

    with ds_store.DSStore.open(filePath, 'w+') as d:
        d['.']['vSrn'] = ('long', 1)
        d['.']['icvl'] = ('type', 'icnv')
        d['.']['bwsp'] = getWindowData_bwsp()
        d['.']['icvp'] = getBackgroundImage_icvp(volume_name)

        for key, iconPosition in iconPositions.items():
            if key == 'app':
                name = app_name
            else:
                name = key
            d[name]['Iloc'] = iconPosition

if __name__ == '__main__':
    if len(sys.argv) != 4:
        print(
            "Usage: %s TARGET_DIR VOLUME_NAME APP_NAME.app" % sys.argv[0],
            file=sys.stderr)
        print("       creates .DS_Store", file=sys.stderr)
        sys.exit(1)
        
    createDSStore(
        target_dir = sys.argv[1],
        volume_name = sys.argv[2],
        app_name = sys.argv[3])
