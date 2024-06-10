# sage_setup: distribution = sagemath-repl
r"""
Sage Commandline Prompts
"""

# ****************************************************************************
#       Copyright (C) 2016 Volker Braun <vbraun.name@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from pygments.token import Token
from IPython.terminal.prompts import Prompts


class SagePrompts(Prompts):

    def in_prompt_tokens(self, cli=None):
        return [
            (Token.Prompt, 'sage: '),
        ]

    def continuation_prompt_tokens(self, cli=None, width=None):
        return [
            (Token.Prompt, '....: '),
        ]

    def rewrite_prompt_tokens(self):
        return [
            (Token.Prompt, '----> '),
        ]

    def out_prompt_tokens(self):
        return [
            (Token.OutPrompt, ''),
        ]


class InterfacePrompts(Prompts):

    def __init__(self, interface_name):
        self.__name = interface_name
        self.__width = len(interface_name)

    def in_prompt_tokens(self, cli=None):
        return [
            (Token.Prompt, self.__name + ': '),
        ]

    def continuation_prompt_tokens(self, cli=None, width=None):
        return [
            (Token.Prompt, '.' * self.__width + ': '),
        ]

    def rewrite_prompt_tokens(self):
        return [
            (Token.Prompt, '-' * self.__width + '> '),
        ]

    def out_prompt_tokens(self):
        return [
            (Token.OutPrompt, ''),
        ]


class DebugPrompts(Prompts):

    def in_prompt_tokens(self, cli=None):
        return [
            (Token.Prompt, 'debug: '),
        ]

    def continuation_prompt_tokens(self, cli=None, width=None):
        return [
            (Token.Prompt, '.....: '),
        ]

    def rewrite_prompt_tokens(self):
        return [
            (Token.Prompt, '-----> '),
        ]

    def out_prompt_tokens(self):
        return [
            (Token.OutPrompt, ''),
        ]
