# Makefile for building the Mac App
#
# This assumes that you have already built Sage. This Makefile then
# compiles the Mac App UI and combines it with the existing Sage
# install into a DMG file.
#
# The resulting DMG is as relocatable as the original Sage install,
# that is, not at all by default. It will contain paths pointing to
# the Sage tree that you started with (OK if you are only debugging
# the Mac App, just don't delete the original Sage directory tree).
#
# To build a relocatable Sage binary you have to use
# https://github.com/sagemath/binary-pkg


ifdef SAGE_ROOT
    $(error This makefile must be run outside of a sage shell)
endif

SAGE_ROOT:=$(shell cd ../.. && pwd)

ifndef SAGE_VERSION
    SAGE_VERSION:=$(shell source ../bin/sage-version.sh && echo $$SAGE_VERSION)
endif

# Meaningful values, on Intel, are "i386" and "x86_64". To prepare a
# 64-bit binary distribution on an older 64-bit OSX machine that boots
# into a 32-bit system, one would set SAGE_APP_TARGET_ARCH=x86_64
ifndef SAGE_APP_TARGET_ARCH
    SAGE_APP_TARGET_ARCH:=$(shell uname -m)
endif

OSX_VERSION:=$(shell sw_vers -productVersion)

TARGET:=sage-$(SAGE_VERSION)-OSX_$(OSX_VERSION)_$(SAGE_APP_TARGET_ARCH)
APP_NAME:=SageMath-$(SAGE_VERSION).app
APP:=$(TARGET)/$(APP_NAME)
DMG:=$(TARGET).app.dmg

CONFIGURATION=Debug
SDKROOT=

SOURCE=*.h *.m

default: $(DMG)
app: $(APP)/Contents/Resources/sage

# Display the variables
debug:
	@echo SAGE_ROOT=$(SAGE_ROOT)
	@echo SAGE_VERSION=$(SAGE_VERSION)
	@echo OSX_VERSION=$(OSX_VERSION)
	@echo SAGE_APP_TARGET_ARCH=$(SAGE_APP_TARGET_ARCH)
	@echo CONFIGURATION=$(CONFIGURATION)
	@echo APP=$(APP)

# Build the Mac App UI binary (not: Sage itself)
$(APP): $(SOURCE)
	mkdir -p $(TARGET)
	xcodebuild -target Sage -configuration $(CONFIGURATION) \
	    ARCHS=$(SAGE_APP_TARGET_ARCH) \
	    ONLY_ACTIVE_ARCH=NO \
	    SDKROOT=$(SDKROOT)
	mkdir -p $(APP)
	cp -pRL build/$(CONFIGURATION)/Sage.app/ $(APP)
	@# Info.plist is a binary plist, so convert it for processing with sed.
	@# I would just change it to be an xml plist, but xcode changes it back.
	plutil -convert xml1 $(APP)/Contents/Info.plist
	sed -i '' "s/SAGE_VERSION/$(SAGE_VERSION)/" $(APP)/Contents/Info.plist

# Copy Sage into the App directory structure
$(APP)/Contents/Resources/sage: $(APP)
	mkdir -p $(APP)/Contents/Resources/sage
	rsync -av --progress \
	    $(SAGE_ROOT)/ $@ \
	    --exclude .git/ \
	    --exclude upstream/ \
	    --exclude src/mac-app/build \
	    --exclude src/mac-app/$(TARGET)
	@# We touch the director here so that we don't keep rsyncing just
	@# because $(APP) is newer than (the rsync preserves time)
	test -d $@ && touch $@

# The background image for the DMG.
# It contains the arrow pointing from the .app to the Applications
# folder.
#
# The background image should be at 72dpi. if a different resolution is
# set in the metadata, Finder produces weird results.
#
$(TARGET)/.background/sage.png: sage-dmg-background.png
	mkdir -p $(TARGET)/.background
	cp sage-dmg-background.png $(TARGET)/.background/sage.png

# Link to Applications
$(TARGET)/Applications:
	rm -f $(TARGET)/Applications
	ln -s /Applications/ $(TARGET)/Applications

# Add the readme file
$(TARGET)/README.txt: $(SAGE_ROOT)/src/bin/sage-README-osx.txt
	cp $^ $@

non_app_files: $(TARGET)/README.txt $(TARGET)/.background/sage.png $(TARGET)/Applications

sage_symlink: $(APP)/Contents/Resources/sage/sage
	rm -f $(APP)/sage
	ln -s Contents/Resources/sage/sage $(APP)/sage

# The arrangement of the icons and the "alias" for the background
# image are stored in a .DS_Store file. The .DS_Store file does
# not store the background image by relative path, but by "alias"
# which contains the volume name as well - so we cannot just
# generate one .DS_Store and copy it from version to version since
# the volume name contains the version.
#
# It would be nice if we could generate the .DS_Store directly
# instead of mounting the image before compressing it.
# Unfortunately, the .DS_Store format is proprietary and changed
# in Mac OS 10.12, so trying this with the ds_store python package
# did not work out (https://bitbucket.org/al45tair/ds_store/issues/7).
# Also see Ticket 20119 and 22739.

# Dependencies need to be in that order because non_app_files
# needs the directory created by $(APP)/Contents/Resources/sage
tmp-$(DMG): $(APP)/Contents/Resources/sage sage_symlink non_app_files
	@# Create initial DMG
	hdiutil create -srcfolder $(TARGET) -volname $(TARGET) -format UDRW tmp-$(DMG)

	@# Mount it so that we can run the Apple Script file on it
	@# to arrange the icons and set the background image
	mkdir -p mnt
	hdiutil attach tmp-$(DMG) -mountpoint mnt
	osascript arrangeIcons.applescript mnt $(APP_NAME)
	hdiutil detach mnt

# Compress the built directory into a DMG
$(DMG): tmp-$(DMG)
	hdiutil convert -format UDBZ -o $(DMG) tmp-$(DMG)	

.PHONY: clean default debug app non_app_files

clean:
	rm -rf $(TARGET) build $(DMG) tmp-$(DMG) mnt
