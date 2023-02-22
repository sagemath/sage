#!/usr/bin/python
# -*- coding: utf-8 -*-

r"""
Python script to sync labels that are migrated from Trac selection lists.
"""

##############################################################################
#       Copyright (C) 2023 Sebastian Oehms <seb.oehms@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
##############################################################################

import os
import sys
from logging import info, warning, getLogger, INFO
from json import loads
from enum import Enum

class SelectionList(Enum):
    """
    Abstract Enum for selection lists.
    """
    def succ(self):
        """
        Return the successor of `self`.
        """
        l = list(self.__class__)
        i = l.index(self) + 1
        if i >= len(l):
            return None
        return l[i]

class ReviewDecision(Enum):
    """
    Enum for `gh pr view` results for `reviewDecision`.
    """
    changes_requested = 'CHANGES_REQUESTED'
    approved = 'APPROVED'
    unclear = 'COMMENTED'

class Priority(SelectionList):
    """
    Enum for priority lables.
    """
    blocker = 'p: blocker /1'
    critical = 'p: critical /2'
    major = 'p: major /3'
    minor = 'p: minor /4'
    trivial = 'p: trivial /5'

class State(SelectionList):
    """
    Enum for state lables.
    """
    positive_review = 's: positive review'
    needs_work = 's: needs work'
    needs_review = 's: needs review'
    needs_info = 's: needs info'


def selection_list(label):
    """
    Return the selection list to which `label` belongs to.
    """
    for sel_list in [Priority, State]:
        for item in sel_list:
            if label == item.value:
                return sel_list
    return None

class GhLabelSynchronizer:
    """
    Handler for access to GitHub issue via the `gh` in the bash command line
    of the GitHub runner.
    """
    def __init__(self, url, actor):
        """
        Python constructor sets the issue / PR url and list of active labels.
        """
        self._url = url
        self._actor = actor
        self._labels = None
        self._author = None
        self._draft = None
        self._open = None
        self._review_decision = None
        self._reviews = None
        self._commits = None

        number = os.path.basename(url)
        self._pr = True
        self._issue = 'pull request #%s' % number
        if url.rfind('issue') != -1:
            self._issue = 'issue #%s' % number
            self._pr = False
        info('Create label handler for %s and actor %s' % (self._issue, self._actor))

    def is_pull_request(self):
        """
        Return if we are treating a pull request.
        """
        return self._pr

    def view(self, key):
        """
        Return data obtained from `gh` command `view`.
        """
        issue = 'issue'
        if self._pr:
            issue = 'pr'
        cmd = 'gh %s view %s --json %s' % (issue, self._url, key)
        from subprocess import check_output
        return loads(check_output(cmd, shell=True))[key]

    def is_open(self):
        """
        Return if the issue res. PR is open.
        """
        if self._open is not None:
            return self._open
        if self.view('state') == 'OPEN':
            self._open = True
        else:
            self._open = False
        info('Issue %s is open %s' % (self._issue, self._open))
        return self._open

    def is_draft(self):
        """
        Return if the PR is a draft.
        """
        if self._draft is not None:
            return self._draft
        if self.is_pull_request():
            self._draft = self.view('isDraft')
        else:
            self._draft = False
        info('Issue %s is draft %s' % (self._issue, self._draft))
        return self._draft

    def get_labels(self):
        """
        Return the list of labels of the issue resp. PR.
        """
        if self._labels is not None:
            return self._labels
        data = self.view('labels')
        self._labels = [l['name'] for l in data]
        info('List of labels for %s: %s' % (self._issue, self._labels))
        return self._labels

    def get_author(self):
        """
        Return the author of the issue resp. PR.
        """
        if self._author is not None:
            return self._author
        data = self.view('author')
        self._author = self.view('author')['login']
        info('Author of %s: %s' % (self._issue, self._author))
        return self._author

    def get_review_decision(self):
        """
        Return the reviewDecision of the PR.
        """
        if not self.is_pull_request():
            return None

        if self._review_decision is not None:
            return self._review_decision

        data = self.view('reviewDecision')
        if data:
            self._review_decision = ReviewDecision(data)
        else:
            self._review_decision = ReviewDecision.unclear
        info('Review decision for %s: %s' % (self._issue, self._review_decision.value))
        return self._review_decision

    def get_reviews(self):
        """
        Return the list of reviews of the PR.
        """
        if not self.is_pull_request():
            return None

        if self._reviews is not None:
            return self._reviews

        self._reviews = self.view('reviews')
        info('Reviews for %s: %s' % (self._issue, self._reviews))
        return self._reviews

    def get_commits(self):
        """
        Return the list of commits of the PR.
        """
        if not self.is_pull_request():
            return None

        if self._commits is not None:
            return self._commits

        self._commits = self.view('commits')
        info('Commits for %s: %s' % (self._issue, self._commits))
        return self._commits

    def gh_cmd(self, cmd, arg, option):
        """
        Perform a system call to `gh` for `cmd` to an isuue resp. PR.
        """
        issue = 'issue'
        if self._pr:
            issue = 'pr'
        cmd_str = 'gh %s %s %s %s "%s"' % (issue, cmd, self._url, option, arg)
        os.system(cmd_str)

    def edit(self, arg, option):
        """
        Perform a system call to `gh` to edit an issue resp. PR.
        """
        self.gh_cmd('edit', arg, option)

    def review(self, arg, text):
        """
        Perform a system call to `gh` to review a PR.
        """
        self.gh_cmd('review', arg, '-b %s' % text)

    def active_partners(self, item):
        """
        Return the list of other labels from the selection list
        of the given one that are already present on the issue / PR.
        """
        sel_list = type(item)
        partners = [i for i in sel_list if i != item and i.value in self.get_labels()]
        info('Active partners of  %s: %s' % (item, partners))
        return partners

    def needs_work(self, only_actor=True):
        """
        Return `True` if the PR needs work. This is the case if
        the review decision requests changes or if there is any
        review reqesting changes.
        """
        ch_req = ReviewDecision.changes_requested
        rev_dec =  self.get_review_decision()
        if rev_dec:
            if rev_dec == ch_req:
               info('PR %s needs work (by decision)' % self._issue)
               return True
            else:
               info('PR %s doesn\'t need work (by decision)' % self._issue)
               return False

        revs = self.get_reviews()
        if only_actor:
            revs = [rev for rev in revs if rev['author']['login'] == self._actor]
        if any(rev['state'] == ch_req.value for rev in revs):
            info('PR %s needs work' % self._issue)
            return True
        info('PR %s doesn\'t need work' % self._issue)
        return False

    def positive_review(self, only_actor=True):
        """
        Return `True` if the PR has positive review. This is the
        case if the review decision is approved or if there is any
        approved review but no changes requesting one.
        """
        appr = ReviewDecision.approved
        rev_dec =  self.get_review_decision()
        if rev_dec:
            if rev_dec == appr:
                info('PR %s has positve review (by decision)' % self._issue)
                return True
            else:
                info('PR %s doesn\'t have positve review (by decision)' % self._issue)
                return False

        if self.needs_work():
            info('PR %s doesn\'t have positve review (needs work)' % self._issue)
            return False

        revs = self.get_reviews()
        if only_actor:
            revs = [rev for rev in revs if rev['author']['login'] == self._actor]
            if any(rev['state'] == appr.value for rev in revs):
                info('PR %s has positve review (by decision)' % self._issue)
                return True
        info('PR %s doesn\'t have positve review (needs work)' % self._issue)
        return False

    def approve_allowed(self):
        """
        Return if the actor has permission to approve this PR.
        """
        author = self.get_author()
        revs = self.get_reviews()
        if not any(rev['author']['authorAssociation'] == 'MEMBER' for rev in revs):
            info('PR %s can\'t be approved because of missing member review' % (self._issue))
            return False

        revs = [rev for rev in revs if rev['author']['login'] != self._actor]
        ch_req = ReviewDecision.changes_requested
        if any(rev['state'] == ch_req.value for rev in revs):
            info('PR %s can\'t be approved by %s since others reqest changes' % (self._issue, self._actor))
            return False

        if author != self._actor:
            info('PR %s can be approved by %s' % (self._issue, self._actor))
            return True

        revs = [rev for rev in revs if rev['author']['login'] != 'github-actions']
        if not revs:
            info('PR %s can\'t be approved by the author %s since no other person reviewed it' % (self._issue, self._actor))
            return False

        comts = self.get_commits()
        authors = sum(com['authors'] for com in comts)
        authors = [auth for auth in authors if not auth['login'] in (self._actor, 'github-actions')]
        if not authors:
            info('PR %s can\'t be approved by the author %s since no other person commited to it' % (self._issue, self._actor))
            return False

        info('PR %s can be approved by the author %s as co-author' % (self._issue, self._actor))
        return True

    def approve(self):
        """
        Approve the PR by the actor.
        """
        self.review('--approve', '%s approved this PR' % self._actor)
        info('PR %s approved by %s' % (self._issue, self._actor))

    def request_changes(self):
        """
        Request changes for this PR by the actor.
        """
        self.review('--request-changes', '%s requested changes for this PR' % self._actor)
        info('Changes requested for PR %s by %s' % (self._issue, self._actor))

    def add_comment(self, text):
        """
        Perform a system call to `gh` to add a comment to an issue or PR.
        """

        self.gh_cmd('comment', text, '-b')
        info('Add comment to %s: %s' % (self._issue, text))

    def add_label(self, label):
        """
        Add the given label to the issue or PR.
        """
        if not label in self.get_labels():
            self.edit(label, '--add-label')
            info('Add label to %s: %s' % (self._issue, label))

    def add_default_label(self, item):
        """
        Add the given label if there is no active partner.
        """
        if not self.active_partners(item):
            self.add_label(item.value)

    def on_label_add(self, label):
        """
        Check if the given label belongs to a selection list. If so, remove
        all other labels of that list. In case of a state label reviews are
        booked accordingly.
        """
        sel_list = selection_list(label)
        if not sel_list:
            return

        item = sel_list(label)
        if label not in self.get_labels():
            # this is possible if two labels of the same selection list
            # have been added in one step (via multiple selection in the
            # pull down menue). In this case `label` has been removed
            # on the `on_label_add` of the first of the two labels
            partn = self.active_partners(item)
            if partn:
                self.add_comment('Label *%s* can not be added due to *%s*!' % (label, partn[0].value))
            else:
                warning('Label %s of %s not found!' % (label, self._issue))
            return

        if sel_list is State:
            if not self.is_pull_request():
                if item != State.needs_info:
                    self.add_comment('Label *%s* can not be added to an issue. Please use it on the corresponding PR' % label)
                    self.remove_label(label)
                    return

            if item == State.positive_review:
                if self.approve_allowed():
                    self.approve()
                elif self.needs_work():
                    # PR still needs work
                    self.add_comment('Label *%s* can not be added. Please use the corresponding functionality of GitHub' % label)
                    self.select_label(State.needs_work)
                    return
                else:
                    # PR still needs review
                    self.add_comment('Label *%s* can not be added. Please use the corresponding functionality of GitHub' % label)
                    if self.is_draft():
                        self.remove_label(label)
                    else:
                        self.select_label(State.needs_review)
                    return

            if item == State.needs_work:
                self.request_changes()
                if not self.needs_work():
                    # change request of actor could not be set
                    self.add_comment('Label *%s* can not be added. Please use the corresponding functionality of GitHub' % label)
                    if self.is_draft():
                        self.remove_label(label)
                    else:
                        self.select_label(State.needs_review)
                    return

        for other in sel_list:
            if other != item:
                self.remove_label(other.value)

    def select_label(self, item):
        """
        Add the given label and remove all others.
        """
        self.add_label(item.value)
        sel_list = type(item)
        for other in sel_list:
            if other != item:
                self.remove_label(other.value)

    def remove_label(self, label):
        """
        Remove the given label from the issue or PR of the handler.
        """
        if label in self.get_labels():
            self.edit(label, '--remove-label')
            info('Remove label from %s: %s' % (self._issue, label))

    def on_label_remove(self, label):
        """
        Check if the given label belongs to a selection list. If so, the
        successor of the label is added except there is none or there
        exists another label of the corresponding list. In case of a
        state label reviews are booked accordingly.
        """
        sel_list = selection_list(label)
        if not sel_list:
            return

        item = sel_list(label)
        if label in self.get_labels():
            # this is possible if two labels of the same selection list
            # have been removed in one step (via multiple selection in the
            # pull down menue). In this case `label` has been added
            # on the `on_label_remove` of the first of the two labels.
            partn = self.active_partners(item)
            if partn:
                self.on_label_add(partn[0].value)
            return

        if sel_list is State:
            if not self.is_pull_request():
                return
            if item == State.positive_review:
                if self.positive_review():
                    self.request_changes()
                    self.select_label(State.needs_work)
                elif self.positive_review(only_actor=False):
                    self.add_comment('Label *%s* can not be removed. Please use the corresponding functionality of GitHub' % label)
                    self.select_label(item)
                elif not self.is_draft():
                    self.select_label(State.needs_review)
                return
            elif item == State.needs_work:
                if self.needs_work(only_actor=False):
                    self.add_comment('Label *%s* can not be removed. Please use the corresponding functionality of GitHub' % label)
                    self.select_label(item)
                elif not self.is_draft():
                    self.select_label(State.needs_review)
                return

        if not self.active_partners(item):
            # if there is no other in the same selection list
            # add the next weaker label if it exists
            succ = sel_list(label).succ()
            if succ:
                self.select_label(succ)
            
    def on_converted_to_draft(self):
        """
        Remove all state labels.
        """
        for item in State:
            self.remove_label(item.value)


###############################################################################
# Main
###############################################################################
cmdline_args = sys.argv[1:]

getLogger().setLevel(INFO)
info('cmdline_args (%s) %s' % (len(cmdline_args), cmdline_args))

if len(cmdline_args) < 4:
    print('Need 5 arguments: action, url, actor, label, rev_state' )
    exit
else:
    action, url, actor, label, rev_state = cmdline_args

info('action: %s' % action)
info('url: %s' % url)
info('label: %s' % label)
info('actor: %s' % actor)
info('rev_state: %s' % rev_state)

gh = GhLabelSynchronizer(url, actor)

if action in ('opened', 'reopened'):
    if gh.is_pull_request():
        if not gh.is_draft():
            gh.add_default_label(State.needs_review)

if action in ('closed', 'reopened'):
    for lab in State:
        gh.remove_label(lab.value)

if action == 'labeled':
    gh.on_label_add(label)

if action == 'unlabeled':
    gh.on_label_remove(label)

if action == 'ready_for_review':
    gh.select_label(State.needs_review)

if action == 'converted_to_draft':
    gh.on_converted_to_draft()

if action == 'submitted':
    if rev_state == 'approved':
        if gh.positive_review():
            gh.select_label(State.positive_review)

    if rev_state == 'changes_requested':
        if gh.needs_work():
            gh.select_label(State.needs_work)

if action in ('review_requested', 'ready_for_review'):
    gh.select_label(State.needs_review)
