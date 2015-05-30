//
//  AppController.m
//  SageMenu
//
//  Created by Ivan Andrus on 19/6/10.
//  Copyright 2010 __MyCompanyName__. All rights reserved.
//

#import "AppController.h"
#import "MyDocument.h"
#import "InputPanelController.h"
#import <WebKit/WebFrame.h>


@implementation AppController

// With help from
// http://www.sonsothunder.com/devres/revolution/tutorials/StatusMenu.html

- (void) awakeFromNib{

    // Used to detect where our files are
    NSBundle *bundle = [NSBundle mainBundle];

    // Allocate and load the images into the application which will be used for our NSStatusItem
    statusImageBlue  = [[NSImage alloc] initWithContentsOfFile:[bundle pathForResource:@"sage-small-blue"  ofType:@"png"]];
    statusImageRed   = [[NSImage alloc] initWithContentsOfFile:[bundle pathForResource:@"sage-small-red"   ofType:@"png"]];
    statusImageGrey  = [[NSImage alloc] initWithContentsOfFile:[bundle pathForResource:@"sage-small-grey"  ofType:@"png"]];
    statusImageGreen = [[NSImage alloc] initWithContentsOfFile:[bundle pathForResource:@"sage-small-green" ofType:@"png"]];

    // Access to the user's defaults
    defaults = [NSUserDefaults standardUserDefaults];

    // Find sageBinary etc.
    [self setupPaths];
    [self ensureReadWrite];

    // Initialize the StatusItem if desired.
    // If we are on Tiger, then showing in the dock doesn't work
    // properly, hence pretend they didn't want it.
    myIsInDock = [defaults boolForKey:@"myShowInDock"] && ![self isTigerOrLess];
    haveStatusItem   = !myIsInDock || [defaults boolForKey:@"alsoShowMenuExtra"];
    useSystemBrowser = !myIsInDock || [defaults boolForKey:@"useSystemBrowser"];
    if ( haveStatusItem ) {
        // Create the NSStatusBar and set its length
        statusItem = [[[NSStatusBar systemStatusBar] statusItemWithLength:NSSquareStatusItemLength] retain];

        // Set the image in NSStatusItem
        [statusItem setImage:statusImageGrey];

        // Tell NSStatusItem what menu to load
        [statusItem setMenu:statusMenu];
        // Set the tooptip for our item
        [statusItem setToolTip:@"Control Sage Notebook Server"];
        // Enable highlighting when menu is opened
        [statusItem setHighlightMode:YES];
    } else {
        [statusItem setEnabled:NO];
    }

    // indicate that we haven't started the server yet
    port = 0;
    neverOpenedFileBrowser = YES;
    URLQueue = [[NSMutableArray arrayWithCapacity:3] retain];

    // Start the sage server, or check if it's running
    if ( [defaults boolForKey:@"startServerOnLaunch"] ) {
        [self startServer:self];
    } else {
        [self serverIsRunning:NO];
    }

    // Set up notifications when an NSTask finishes.
    // For us this will be for checking if the server is running
    [[NSNotificationCenter defaultCenter] addObserver:self
                                             selector:@selector(taskTerminated:)
                                                 name:NSTaskDidTerminateNotification
                                               object:nil];
}

- (void) dealloc {
    // Release the images we loaded into memory
    [statusImageBlue  release];
    [statusImageRed   release];
    [statusImageGrey  release];
    [statusImageGreen release];
    [sageBinary release];
    [logPath release];
    [theTask release];
    [taskPipe release];
    [URLQueue release];
    [super dealloc];
}


-(IBAction)startServer:(id)sender{
    // TODO: Check to see if it's running before attempting to start
    NSLog(@"Starting server");
    if (haveStatusItem)  [statusItem setImage:statusImageGreen];
    NSString *scriptPath  = [[NSBundle mainBundle] pathForResource:@"start-sage" ofType:@"sh"];

    // Add SAGE_BROWSER to environment to point back to this application
    if ( !useSystemBrowser ) {
        NSString *browserPath = [[NSBundle mainBundle] pathForResource:@"open-location" ofType:@"sh"];
        setenv("SAGE_BROWSER", [browserPath UTF8String], 1); // this overwrites, should it?
    }

    // Create a task to start the server
    
    // Get any default options they might have for this session
    [defaults synchronize];
    NSString *defArgs = [[defaults dictionaryForKey:@"DefaultArguments"]
                         objectForKey:@"notebook"];
    launchTask = [[NSTask launchedTaskWithLaunchPath:scriptPath
                                           arguments:[NSArray arrayWithObjects:sageBinary,
                                                      logPath,
                                                      defArgs, // May be nil, but that's okay
                                                      nil]]
                  retain];

    // Open loading page since it can take a while to start
    [self browseRemoteURL:[[NSBundle mainBundle] pathForResource:@"loading-page" ofType:@"html"]];

    // Get info about the server if we're not going to get it via opening a page
    if ( useSystemBrowser ) {
        [self serverIsRunning:YES];
    }
}

-(BOOL)serverIsRunning:(BOOL)wait{

    // Start busy polling until the server starts
    if ( theTask == nil && taskPipe == nil ) {
        theTask  = [[NSTask alloc] init];
        taskPipe = [[NSPipe alloc] init];
        [theTask setStandardOutput:taskPipe];
        [theTask setLaunchPath:[[NSBundle mainBundle] pathForResource:@"sage-is-running-on-port" ofType:@"sh"]];
        if (wait)  [theTask setArguments:[NSArray arrayWithObject:@"--wait"]];
        [theTask launch];
    }
    return NO;
}

-(void)serverStartedWithPort:(int)p{
    if (haveStatusItem)  [statusItem setImage:statusImageBlue];
    port = p;
    if ( [URLQueue count] > 0 ) {
        NSEnumerator *e = [URLQueue objectEnumerator];
        id url;
        while (url = [e nextObject]) {
            [self browseLocalSageURL:url];
        }
        [URLQueue removeAllObjects];
    }
}

- (void)taskTerminated:(NSNotification *)aNotification {

    NSTask *theObject = [aNotification object];
    if (theObject == theTask) {
        const int status = [theObject terminationStatus];
        if (status == 0) {
            // Parse the output
            NSData *data = [[taskPipe fileHandleForReading] readDataToEndOfFile];
            NSString* s = [[NSString alloc] initWithBytes:[data bytes]
                                                   length:[data length]
                                                 encoding:NSUTF8StringEncoding];
            const int p = [s intValue];
            [s release];
            [self serverStartedWithPort:p];
        } else {
            // We failed, so tell the user
            if (haveStatusItem) {
                [statusItem setImage:statusImageGrey];
            }
            port = 0;
        }
        // Reset for next time.
        [theTask release];
        theTask = nil;
        [taskPipe release];
        taskPipe = nil;
    } else if (theObject == launchTask ) {
        
        const int status = [theObject terminationStatus];
        if (status != 0) {
            // We failed, so tell the user
            if (haveStatusItem) {
                [statusItem setImage:statusImageGrey];
            }
            port = 0;
            NSAlert *alert = [NSAlert alertWithMessageText:@"Sage Server failed to start"
                                             defaultButton:@"View Log"
                                           alternateButton:@"Cancel"
                                               otherButton:nil
                                 informativeTextWithFormat:@"For some reason the Sage server failed to start.  "
                              "Please check the log for clues, and have that information handy when asking for help."];
            [alert setAlertStyle:NSWarningAlertStyle];
            NSInteger resp = [alert runModal];
            if (resp == NSAlertDefaultReturn) {
                // View Log
                [self viewSageLog:self];
            } else {
                // Cancel
            }
        }
        // Reset for next time.
        [launchTask release];
        launchTask = nil;
        
    } else {
        // NSLog(@"Got called for a different task.");
    }
}

-(IBAction)stopServer:(id)sender{
    if (haveStatusItem)  [statusItem setImage:statusImageRed];

    // Get the pid of the Sage server
    NSString *pidFile = [@"~/.sage/sage_notebook.sagenb/twistd.pid" stringByStandardizingPath];
    NSString *pid = [NSString stringWithContentsOfFile:pidFile
                                              encoding:NSUTF8StringEncoding
                                                 error:NULL];

    if (pid == nil) {
        // Get the pid of the Sage server
        pidFile = [@"~/.sage/sage_notebook.sagenb/sagenb.pid" stringByStandardizingPath];
        pid = [NSString stringWithContentsOfFile:pidFile
                                        encoding:NSUTF8StringEncoding
                                           error:NULL];
    }

    NSLog(@"Stopping server with pid: %@", pid );
    if (pid != nil) {
        kill([pid intValue], SIGTERM);
    }

    if (haveStatusItem)  [statusItem setImage:statusImageGrey];
    port = 0;
}

// To create an alternate menu, in IB create another menu item, give it a key equivalent of opt/alt and check the alternate box (left most tab of inspector)
-(IBAction)stopServerAndQuit:(id)sender{

    [self stopServer:self];

    // Tell the application to quit
    [NSApp performSelector:@selector(terminate:) withObject:nil afterDelay:0.0];
}

-(IBAction)viewSageLog:(id)sender{
    if (logPath != nil) {
        // open files with the default viewer (I think the default is Console.app)
        // http://lethain.com/entry/2008/apr/05/opening-files-with-associated-app-in-cocoa/
        NSWorkspace *workspace = [NSWorkspace sharedWorkspace];
        [workspace openFile:logPath];
    }
}

-(void)setupPaths{

    // Find the log path
    NSFileManager *fileMgr = [NSFileManager defaultManager];
    NSArray *directories = NSSearchPathForDirectoriesInDomains(NSLibraryDirectory, NSUserDomainMask, YES);
    NSString *tmpLogPath;
    if ([directories count] > 0) {
        tmpLogPath = [directories objectAtIndex:0];
        tmpLogPath = [tmpLogPath stringByAppendingPathComponent:@"Logs"];

        if ( [fileMgr fileExistsAtPath:tmpLogPath] || [fileMgr createDirectoryAtPath:tmpLogPath attributes:nil] ) {
            logPath = [tmpLogPath stringByAppendingPathComponent:@"sage.log"];
/*  If we want to send our log there...
            int fd = open([logPath fileSystemRepresentation], (O_RDWR|O_CREAT|O_TRUNC), (S_IRWXU|S_IRWXG|S_IRWXO));
            if (fd != -1) {
                result = asl_add_log_file(NULL, fd);
            }
 */
        } else {
            logPath = [[NSBundle mainBundle] pathForResource:@"sage" ofType:@"log"];
            NSLog(@"Couldn't create the directory (%@) for log file.  Going to log to %@.",tmpLogPath,logPath);
        }
    } else {
        logPath = [[NSBundle mainBundle] pathForResource:@"sage" ofType:@"log"];
        NSLog(@"Something is fishy: couldn't find a path for log files.  Going to log to %@.", logPath);
    }
    [logPath retain];

    // ### Find the sage binary ###

    // If they have a plist entry telling where it is try that.
    sageBinary = [defaults objectForKey:@"SageBinary"];
    // If that isn't wanted or isn't executable, try a sage built in to the application
    BOOL isDir = YES;
    // If the file is a directory, see if it's SAGE_ROOT
    if ( [fileMgr fileExistsAtPath:sageBinary isDirectory:&isDir] && isDir ) {
        [defaults setObject:[sageBinary stringByAppendingPathComponent:@"sage"]
                     forKey:@"SageBinary"];
        sageBinary = [defaults objectForKey:@"SageBinary"];
    }
    // Put isDir last since technically it's value is undefined if the file doesn't exist
    if ( ![defaults boolForKey:@"useAltSageBinary"] || ![fileMgr isExecutableFileAtPath:sageBinary] ) {
        NSString * path = [[NSBundle mainBundle] pathForResource:@"sage" ofType:nil inDirectory:@"sage"];
        sageBinary = path ? [[NSString alloc] initWithString:path] : nil;
        [defaults setBool:NO forKey:@"useAltSageBinary"];
    }

    // If that doesn't work then have them locate a binary for us
    if ( !sageBinary || ![fileMgr isExecutableFileAtPath:sageBinary] ) {

        // Create a File Open Dialog class
        NSOpenPanel *openDlg = [NSOpenPanel openPanel];

        // Enable the selection of files and directories
        [openDlg setTitle:@"Please choose a Sage executable"];
        [openDlg setMessage:@"This application did not come with a Sage distribution, and there is no valid alternative specified.\n\
Please choose a Sage executable to use from now on.  If you do not, sage is assumed to be in PATH.\n\
You can change it later in Preferences."];
        [openDlg setCanChooseFiles:YES];
        [openDlg setCanChooseDirectories:YES];

        // Display the dialog.  If the OK button was pressed,
        // process the files.
        while ( [openDlg runModalForDirectory:nil file:nil] == NSOKButton ) {
            sageBinary = [[openDlg filenames] objectAtIndex:0];
            // if they give a folder, look for sage inside
            if ( [fileMgr fileExistsAtPath:sageBinary isDirectory:&isDir] && isDir ) {
                sageBinary = [sageBinary stringByAppendingPathComponent:@"sage"];
            }
            // Sanity check for the validity of the Sage Binary
            if ( [fileMgr isExecutableFileAtPath:sageBinary] ) {
                // Save for future sessions
                [defaults setBool:YES forKey:@"useAltSageBinary"];
                [defaults setObject:sageBinary forKey:@"SageBinary"];
                [sageBinary retain];
                return;
            }
            [openDlg setMessage:@"That does not appear to be a valid sage executable.\nPlease choose another, or cancel to assume sage is in PATH."];
        }

        // Quit since there's no point going on.
        // [NSApp performSelector:@selector(terminate:) withObject:nil afterDelay:0.0];

        NSLog(@"WARNING: Could not find a good sage executable, falling back to sage and hoping it's in PATH.");
        sageBinary = @"sage";
    }
}

-(void)ensureReadWrite {
    NSFileManager *filemgr = [NSFileManager defaultManager];
    NSLog(@"Checking if sageBinary (%@) is writeable.",sageBinary);
    if ( ! [filemgr isWritableFileAtPath:sageBinary] ) {
        NSAlert *alert = [NSAlert alertWithMessageText:@"Read-only Sage"
                                         defaultButton:@"Quit"
                                       alternateButton:@"Preferences"
                                           otherButton:@"Continue"
                             informativeTextWithFormat:@"You are attempting to run Sage.app with a read-only copy of Sage "
                          "(most likely due to running it from the disk image).  "
                          "Unfortunately, this is not supported for technical reasons.  \n"
                          "Please drag Sage.app to your hard-drive and run it from there, "
                          "or choose a different executable in Preferences."];
        [alert setAlertStyle:NSWarningAlertStyle];
        NSInteger resp = [alert runModal];
        if (resp == NSAlertDefaultReturn) {// Quit
            NSLog(@"Quitting after a read-only Sage warning.");
            [NSApp terminate:self];
        } else if ( resp == NSAlertAlternateReturn) { // Continue
            NSLog(@"Preferences after a read-only Sage warning.");
            [self showPreferences:self];
        } else {
            NSLog(@"Continuing from read-only Sage warning.");
        }
    }
}

-(IBAction)revealInFinder:(id)sender{
    if ( [[sender title] isEqualToString:@"Reveal in Shell"] ) {
        [self terminalRun:[NSString stringWithFormat:@"cd '%@' && $SHELL",
                           [sageBinary stringByDeletingLastPathComponent]]];
    } else {
        [[NSWorkspace sharedWorkspace] selectFile:[sageBinary stringByDeletingLastPathComponent]
                         inFileViewerRootedAtPath:nil];
    }
}

-(IBAction)openNotebook:(id)sender{
    [self browseLocalSageURL:@""];
}

-(IBAction)newWorksheet:(id)sender{
    [self browseLocalSageURL:@"new_worksheet"];
}

-(IBAction)showPreferences:(id)sender{
    [NSApp activateIgnoringOtherApps:YES];
    [prefWindow makeKeyAndOrderFront:self];
}

-(IBAction)browseLocalSageURL:(id)sender{
    NSString *sageURL;
    if ([sender isKindOfClass:[NSString class]]) {
        sageURL = sender;
    } else {
        sageURL = [[defaults arrayForKey:@"sageURLs"] objectAtIndex:[sender tag]];
    }
    // The server is not running
    if ( port == 0 && [defaults boolForKey:@"autoStartServer"] ) {
        // Queue the URL up for opening and start the server
        // Do I need to retain it??
        [URLQueue addObject:sageURL];
        [self startServer:self];
    } else {
        // Browse to the url right away
        [self sageBrowse:[NSString stringWithFormat:@"http://localhost:%d/%@", port, sageURL]];
    }
}

-(IBAction)browseRemoteURL:(id)sender{
    NSString *sageURL;
    if ([sender isKindOfClass:[NSString class]]) {
        sageURL = sender;
    } else {
        sageURL = [[defaults arrayForKey:@"sageURLs"] objectAtIndex:[sender tag]];
    }
    [self sageBrowse:sageURL];
}

-(void)sageBrowse:(NSString*)location{

    if ( !useSystemBrowser ) {
        [[NSApplication sharedApplication] activateIgnoringOtherApps:TRUE];

        NSError *outError = nil;
        id myDocument = [[NSDocumentController sharedDocumentController]
                         openUntitledDocumentAndDisplay:YES error:&outError];
        if ( myDocument == nil ) {
            [NSApp presentError:outError];
            NSLog(@"sageBrowser: Error creating document: %@", [outError localizedDescription]);
        } else {
            [[[myDocument webView] mainFrame]
            loadRequest:[NSURLRequest requestWithURL:[NSURL URLWithString:location]]];
        }

    } else if ( [defaults boolForKey:@"respectSAGE_BROWSER"] ) {

        // TODO: escape quotes in location
        NSString *command = [NSString
                             stringWithFormat:@"%@ -min -c 'import sage.misc.viewer as b; os.system(b.browser() + \" %@\")' &",
                             sageBinary,
                             location];

        // TODO: Should probably make this use NSTask
        system([command UTF8String]);
    } else {

        if ( [location characterAtIndex:0] == '/' ) {
            [[NSWorkspace sharedWorkspace] openFile:location];
        } else {
            [[NSWorkspace sharedWorkspace] openURL:[NSURL URLWithString:location]];
        }
    }
}


-(NSString*)convertMenuTitleToSageCommand:(NSString*)title{

    if ( [title isEqualToString:@"Sage"] || [title isEqualToString:@"Sage (advanced)"] || [title isEqualToString:@"Terminal Session"] ) {
        // A few special cases to open sage itself
        return nil;
    } else if ( ([title length] > 2) && [[NSCharacterSet uppercaseLetterCharacterSet] characterIsMember:[title characterAtIndex:0]] ) {
        // If it's capitalized, and more than one character then use the lowercased first letter.
        // This is so things like Build and Test can work, but R and M2 will still work.
        // This is really a hack, because I'm too lazy to create a bunch of different methods (and I think it's ugly)
        unichar first = [[title lowercaseString] characterAtIndex:0];
        return [NSString stringWithCharacters:&first length:1];
    } else {
        // If it's lowercased, assume it's the command, but remove ... from the end
        return [title stringByTrimmingCharactersInSet:
                [NSCharacterSet characterSetWithCharactersInString:
                 [NSString stringWithFormat:@"%C", ((unsigned short)0x2026)]]]; // @"…"
    }
}

-(IBAction)terminalSession:(id)sender{
    [self sageTerminalRun: [self convertMenuTitleToSageCommand:[sender title]] withArguments: nil];
}

-(IBAction)terminalSessionPromptForInput:(id)sender{

    NSString *sessionType = [self convertMenuTitleToSageCommand:[sender title]];
    NSString *command;
    if ( [sessionType length] > 1 ) {
        command = [sageBinary stringByAppendingFormat:@" --%@", sessionType];
    } else if ( [sessionType length] > 0 ) {
        command = [sageBinary stringByAppendingFormat:@" -%@", sessionType];
    } else {
        command = sageBinary;
    }

    [defaults synchronize];
    NSString *defArgs = [[defaults dictionaryForKey:@"DefaultArguments"]
                         objectForKey:command];

    [inputPanelController runCommand:command
                          withPrompt:[self createPrompt:sessionType forCommand:command]
                       withArguments:defArgs
                      editingCommand:[defaults boolForKey:@"editFullCommands"]];
}

-(NSString*)createPrompt:(NSString*)sessionType forCommand:(NSString*)command{
    return [NSString stringWithFormat:@"Going to run sage %@\nPlease enter any arguments, escaped as you would for a shell.\n\nThe command will be run as\n%@ %C",
            sessionType ? sessionType : @"", command, ((unsigned short)0x2026)];
}

-(IBAction)terminalSessionPromptForFile:(id)sender{

    // Create a File Open Dialog class
    NSOpenPanel *openDlg = [NSOpenPanel openPanel];

    // Enable the selection of files and directories
    [openDlg setCanChooseFiles:YES];
    [openDlg setCanChooseDirectories:YES];
    [openDlg setAllowsMultipleSelection:YES];
    [openDlg setTitle:[NSString stringWithFormat:@"Choose file(s) for %@",[sender title]]];

    // Display the dialog.  If the OK button was pressed,
    // process the files.
    NSString * base_dir = nil;
    if (neverOpenedFileBrowser) {
        base_dir = [NSString stringWithFormat:@"%@/../devel/sage/sage",sageBinary];
        neverOpenedFileBrowser=NO;
    }
    // If they supply files, then run the command
    if ( [openDlg runModalForDirectory:base_dir file:nil] == NSOKButton ) {
        [self sageTerminalRun:[self convertMenuTitleToSageCommand:[sender title]]
                withArguments:[openDlg filenames]];
    }
}

-(void)sageTerminalRun:(NSString*)sessionType withArguments:(NSArray*)arguments{
    NSString *command;
    if ( sessionType == nil ) {
        NSLog(@"starting sage" );
        command = sageBinary;
    } else if ( [sessionType length] > 1 ) {
        command = [sageBinary stringByAppendingFormat:@" --%@", sessionType];
    } else {
        command = [sageBinary stringByAppendingFormat:@" -%@", sessionType];
    }

    // Get any default options they might have for this session
    [defaults synchronize];
    NSString *defArgs = [[defaults dictionaryForKey:@"DefaultArguments"]
                         objectForKey:(sessionType != nil) ? sessionType : @"sage" ];
    if ( defArgs != nil ) {
        command = [command stringByAppendingFormat:@" %@", defArgs];
    }
    if ( arguments != nil ) {
        for( int i = 0; i < [arguments count]; i++ ) {
            command = [command stringByAppendingFormat:@" %@", [arguments objectAtIndex:i]];
        }
    }

    // Hold command key to edit before running
    if ( [defaults boolForKey:@"alwaysPromptForArguments"] || [[NSApp currentEvent] modifierFlags] & NSCommandKeyMask ) {
        [inputPanelController runCommand:command
                              withPrompt:[self createPrompt:sessionType forCommand:command]
                           withArguments:defArgs
                          editingCommand:YES];
    } else {
        [self terminalRun:command];
    }
}

-(void)terminalRun:(NSString*)command{
    NSLog(@"Running command: %@", command);

    // Escape quotes and backslashes in the command
    // I think that's all we need to handle for applescript itself
    NSMutableString * escapedCommand = [NSMutableString stringWithString:command];
    [escapedCommand replaceOccurrencesOfString:@"\\"
                                    withString:@"\\\\"
                                       options:0
                                         range:NSMakeRange(0, [escapedCommand length])];
    [escapedCommand replaceOccurrencesOfString:@"\""
                                    withString:@"\\\""
                                       options:0
                                         range:NSMakeRange(0, [escapedCommand length])];
    // We can't use the (arguably easier) stringByReplacingOccurrencesOfString:withString since that's 10.5+
    // NSString *escapedCommand = [[command stringByReplacingOccurrencesOfString:@"\\" withString:@"\\\\"]
    //                                      stringByReplacingOccurrencesOfString:@"\"" withString:@"\\\""];

    // Which applescript to run
    NSString *ApplescriptKey = [defaults objectForKey:@"TerminalEmulator"];
    // Print the command into the applescript
    NSString *bringAppToFrontScript =
        [NSString stringWithFormat:[[defaults dictionaryForKey:@"TerminalEmulatorList"]
                                    objectForKey:ApplescriptKey],
         escapedCommand];

    // NSLog(@"Executing applescript: %@", bringAppToFrontScript);

    NSDictionary* errorDict;
    NSAppleEventDescriptor *returnDescriptor = NULL;
    NSAppleScript* scriptObject = [[NSAppleScript alloc]
                                   initWithSource:bringAppToFrontScript];
    returnDescriptor = [scriptObject executeAndReturnError: &errorDict];
    if ( returnDescriptor == nil ) {
        NSLog(@"terminalRun: Error running Applescript: %@", errorDict);
    }
    [scriptObject release];
}

// http://www.cocoadev.com/index.pl?DeterminingOSVersion
-(BOOL)isTigerOrLess{
    OSErr err;
    SInt32 version;
    if ((err = Gestalt(gestaltSystemVersionMajor, &version)) != noErr) {
        NSLog(@"Unable to determine gestaltSystemVersionMajor: %hd",err);
        return YES;
    }
    if ( version < 10 ) return YES; // Of course this should never happen...
    if ((err = Gestalt(gestaltSystemVersionMinor, &version)) != noErr) {
        NSLog(@"Unable to determine gestaltSystemVersionMinor: %hd",err);
        return YES;
    }
    if ( version < 5 )  return YES;
    return NO;
}

// TODO: make installing packages easy -- stringByLaunchingPath:withArguments:error:
// TODO: maybe this should be written in py-objc so that we can call into sage directly (but then we would have to worry about environment etc.)
// TODO: make some services (search for NSSendTypes) -- pack/unpack spkg, extract sws from pdf, crap/fixdoctests/preparse/Test/coverage/pkg/pkg_nc/etc.

// TODO: open files such as .sws, .sage, .py, .spkg, -- .pdf (and extract sws from them), .htm, whatever else I can handle
// TODO: quicklook generator, spotlight importer -- use UTI
// NOTE: http://developer.apple.com/mac/library/documentation/Miscellaneous/Reference/UTIRef/Articles/System-DeclaredUniformTypeIdentifiers.html
// TODO: icons for files -- they need some help with the alpha channel.  I clearly don't know what I'm doing.  I should really make them all from by script...


@end
