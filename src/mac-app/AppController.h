//
//  AppController.h
//  SageMenu
//
//  Created by Ivan Andrus on 19/6/10.
//  Copyright 2010 __MyCompanyName__. All rights reserved.
//

#import <Cocoa/Cocoa.h>


@interface AppController : NSObject {

    /* Our outlets which allow us to access the interface */
    IBOutlet NSMenu *statusMenu;
    IBOutlet id appDelegate;
    IBOutlet id prefWindow;
    IBOutlet id inputPanelController;

    NSStatusItem *statusItem;
    NSImage *statusImageBlue;
    NSImage *statusImageGrey;
    NSImage *statusImageGreen;
    NSImage *statusImageRed;

    NSString *sageBinary;
    NSString *logPath;
    NSString *jupyterURL;
    NSMutableArray *URLQueue;

    NSUserDefaults *defaults;

    NSTask *theTask;
    NSTask *launchTask;
    NSPipe *taskPipe;
    NSTask *jupyterTask;

    int port;
    BOOL myIsInDock, haveStatusItem, useSystemBrowser, neverOpenedFileBrowser;

}

// Server control
-(IBAction)startJupyter:(id)sender;
-(IBAction)stopJupyter:(id)sender;
-(void)receivedData:(NSNotification *)notif;
-(IBAction)startServer:(id)sender;
-(IBAction)stopServer:(id)sender;
-(BOOL)serverIsRunning:(BOOL)wait;
-(void)serverStartedWithPort:(int)port;
-(void)taskTerminated:(NSNotification *)aNotification;

// Notebook functions
-(IBAction)openNotebook:(id)sender;
-(IBAction)newWorksheet:(id)sender;
-(IBAction)browseLocalSageURL:(id)sender;
-(IBAction)browseRemoteURL:(id)sender;

// Terminal and advanced
-(IBAction)terminalSession:(id)sender;
-(IBAction)viewSageLog:(id)sender;
-(IBAction)revealInFinder:(id)sender;
-(IBAction)terminalSessionPromptForFile:(id)sender;
-(IBAction)terminalSessionPromptForInput:(id)sender;
-(NSString*)convertMenuTitleToSageCommand:(NSString*)title;

-(IBAction)showPreferences:(id)sender;

-(void)ensureReadWrite;
-(void)offerNotebookUpgrade;
-(IBAction)upgradeNotebook:(id)sender;
-(void)setupPaths;

// Quit
-(IBAction)stopServerAndQuit:(id)sender;

// Ancillary functions
-(void)sageBrowse:(NSString*)location;
-(void)sageTerminalRun:(NSString*)sessionType withArguments:(NSArray*)arguments;
-(void)terminalRun:(NSString*)command;
-(NSString*)createPrompt:(NSString*)sessionType forCommand:(NSString*)command;
-(BOOL)isTigerOrLess;

@end
