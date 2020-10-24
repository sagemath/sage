//
//  MyDocument.m
//  Sage
//
//  Created by Ivan Andrus on 26/6/10.
//  Copyright 2010 __MyCompanyName__. All rights reserved.
//

#import "MyDocument.h"
#import <WebKit/WebFrame.h>
#import <WebKit/WebUIDelegate.h>
#import <WebKit/WebDataSource.h>

@implementation MyDocument

- (id)init{
    self = [super init];
    if (self) {

        // Add your subclass-specific initialization here.
        // If an error occurs here, send a [self release] message and return nil.
        NSNotificationCenter *nc = [NSNotificationCenter defaultCenter];
        [nc addObserver:self selector:@selector(webViewProgressStarted:)  name:WebViewProgressStartedNotification  object:webView];
        [nc addObserver:self selector:@selector(webViewProgressFinished:) name:WebViewProgressFinishedNotification object:webView];

        // We don't want undo
        [self setHasUndoManager:NO];
    }
    return self;
}

- (NSString *)windowNibName{
    // Override returning the nib file name of the document
    // If you need to use a subclass of NSWindowController or if your document supports multiple NSWindowControllers, you should remove this method and override -makeWindowControllers instead.
    return @"MyDocument";
}

- (void)windowControllerDidLoadNib:(NSWindowController *) aController
{
    // TODO: this is slightly underhanded but easier than building our own from scratch.
    [webView setApplicationNameForUserAgent:@"Safari/528.16 SageBrowser"];

    [super windowControllerDidLoadNib:aController];

    // Set up some delegates.  Perhaps this could/should be done in the nib file
    [webView setGroupName:@"MyDocument"];
    [webView setUIDelegate:self];
    [webView setFrameLoadDelegate:self];
}

// TODO: this will allow saving
- (NSData *)dataOfType:(NSString *)typeName error:(NSError **)outError
{
    // Insert code here to write your document to data of the specified type. If the given outError != NULL, ensure that you set *outError when returning nil.

    // You can also choose to override -fileWrapperOfType:error:, -writeToURL:ofType:error:, or -writeToURL:ofType:forSaveOperation:originalContentsURL:error: instead.

    // For applications targeted for Panther or earlier systems, you should use the deprecated API -dataRepresentationOfType:. In this case you can also choose to override -fileWrapperRepresentationOfType: or -writeToFile:ofType: instead.
    NSLog(@"well at least I made it there");

    if ( outError != NULL ) {
        *outError = [NSError errorWithDomain:NSOSStatusErrorDomain code:unimpErr userInfo:NULL];
    }

    return nil;
}


- (BOOL)readFromURL:(NSURL *)absoluteURL ofType:(NSString *)typeName error:(NSError **)outError {

    NSLog(@"well at least I made it to open a url");
    if ( outError != NULL ) {
        *outError = [NSError errorWithDomain:NSOSStatusErrorDomain code:unimpErr userInfo:NULL];
    }
    return YES;
}

- (BOOL)readFromData:(NSData *)data ofType:(NSString *)typeName error:(NSError **)outError
{
    // Insert code here to read your document from the given data of the specified type.  If the given outError != NULL, ensure that you set *outError when returning NO.

    // You can also choose to override -readFromFileWrapper:ofType:error: or -readFromURL:ofType:error: instead.

    // For applications targeted for Panther or earlier systems, you should use the deprecated API -loadDataRepresentation:ofType. In this case you can also choose to override -readFromFile:ofType: or -loadFileWrapperRepresentation:ofType: instead.
    NSLog(@"well at least I made it here");
    if ( outError != NULL ) {
        *outError = [NSError errorWithDomain:NSOSStatusErrorDomain code:unimpErr userInfo:NULL];
    }
    return YES;
}

// From Fluidium
- (BOOL)writeToURL:(NSURL *)absoluteURL ofType:(NSString *)typeName error:(NSError **)outError {
    NSData *archiveData = [[[[webView mainFrame] dataSource] webArchive] data];
    return [archiveData writeToURL:absoluteURL options:0 error:outError];
}

- (id)webView{
    return webView;
}

- (IBAction)connectURL:(id)sender{
    [urlString setStringValue:[sender stringValue]];
    [[webView mainFrame] loadRequest:
     [NSURLRequest requestWithURL:
      [NSURL URLWithString:
       [sender stringValue]]]];
}

- (WebView *)webView:(WebView *)sender createWebViewWithRequest:(NSURLRequest *)request{

    NSError *outError = nil;
    id myDocument = [[NSDocumentController sharedDocumentController]
                     openUntitledDocumentAndDisplay:YES error:&outError];
    if ( myDocument == nil ) {
        [NSApp presentError:outError];
        NSLog(@"sageBrowser: Error creating document: %@", [outError localizedDescription]);
    } else {
        [[[myDocument webView] mainFrame]
         loadRequest:request];

    }

    return [myDocument webView];
}

- (void)webViewShow:(WebView *)sender{
    id myDocument = [[NSDocumentController sharedDocumentController] documentForWindow:[sender window]];
    [myDocument showWindows];
}

- (IBAction)browseURL:(NSString*)theURL{

    id myDocument = [[NSDocumentController sharedDocumentController]
                     openUntitledDocumentOfType:@"DocumentType"
                     display:YES];

    [[[myDocument webView] mainFrame]
     loadRequest:[NSURLRequest requestWithURL:[NSURL URLWithString:theURL]]];
}

- (void)webView:(WebView *)sender didStartProvisionalLoadForFrame:(WebFrame *)frame{
    // Only report feedback for the main frame.
    if (frame == [sender mainFrame]){
        NSString *url = [[[[frame provisionalDataSource] request] URL] absoluteString];
        [urlString setStringValue:url];
    }
}

// Taken from Fluidium

- (void)webView:(WebView *)wv runOpenPanelForFileButtonWithResultListener:(id <WebOpenPanelResultListener>)listener {
    NSOpenPanel *openPanel = [NSOpenPanel openPanel];
    [openPanel beginSheetForDirectory:nil
                                 file:@""
                       modalForWindow:[[self webView] window]
                        modalDelegate:self
                       didEndSelector:@selector(openPanelDidEnd:returnCode:contextInfo:)
                          contextInfo:[listener retain]]; // retained
}

- (void)openPanelDidEnd:(NSSavePanel *)openPanel returnCode:(int)code contextInfo:(id <WebOpenPanelResultListener>)listener {
    [listener autorelease]; // released

    if (NSOKButton == code) {
        [listener chooseFilename:[openPanel filename]];
    }
}

- (IBAction)printDocument:(id)sender{
    [[[[webView mainFrame] frameView] documentView] print:sender];
}

#pragma mark WebProgressNotifications

- (void)webViewProgressStarted:(NSNotification *)n {
//    NSLog(@"progress started: %@", progressIndicator);
    [progressIndicator startAnimation:self];
}


- (void)webViewProgressFinished:(NSNotification *)n {
//    NSLog(@"progress stopped");
    [progressIndicator stopAnimation:self];
}

@end
