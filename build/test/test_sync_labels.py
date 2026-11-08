"""Offline regression tests for the privileged label synchronization script."""

import importlib.util
import json
import logging
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch


source = Path(__file__).resolve().parents[2] / '.github/sync_labels.py'
spec = importlib.util.spec_from_file_location('sync_labels', source)
sync_labels = importlib.util.module_from_spec(spec)
spec.loader.exec_module(sync_labels)


class LabelCommandsTestCase(unittest.TestCase):
    def setUp(self):
        self.addCleanup(logging.getLogger().setLevel, logging.getLogger().level)
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.root = Path(self.directory.name)
        self.calls = self.root / 'calls.json'
        gh = self.root / 'gh'
        gh.write_text(f'''#!{sys.executable}
import json
import os
import sys
from pathlib import Path
Path(os.environ['FAKE_GH_CALLS']).write_text(json.dumps(sys.argv[1:]))
if sys.argv[1] == 'api' or 'view' in sys.argv:
    print(os.environ.get('FAKE_GH_RESPONSE', '{{}}'))
sys.exit(int(os.environ.get('FAKE_GH_EXIT', '0')))
''')
        gh.chmod(0o755)
        self.environment = patch.dict(os.environ, {
            'PATH': str(self.root) + os.pathsep + os.environ.get('PATH', ''),
            'FAKE_GH_CALLS': str(self.calls),
            'FAKE_GH_RESPONSE': '{}',
            'FAKE_GH_EXIT': '0',
        })
        self.environment.start()
        self.addCleanup(self.environment.stop)
        # Avoid the constructor's API queries and comment cleanup.
        self.handler = object.__new__(sync_labels.GhLabelSynchronizer)
        self.handler._url = 'https://github.com/example/project/pull/42'
        self.handler._pr = True
        self.handler._issue = 'pull request #42'
        self.marker = self.root / 'unexpected-shell-execution'
        self.payload = f'Quotes " and \' and $(touch {self.marker}) `touch {self.marker}`\nnext line'

    def arguments(self):
        self.assertFalse(self.marker.exists(), 'Text was evaluated by a shell')
        return json.loads(self.calls.read_text())

    def test_comment_text_is_literal(self):
        self.handler.add_comment(self.payload)
        self.assertEqual(self.arguments(), [
            '--repo', 'https://github.com/example/project', 'pr', 'comment',
            self.handler._url, '-b', self.payload,
        ])

    def test_review_text_is_literal(self):
        self.handler.review('--comment', self.payload)
        self.assertEqual(self.arguments()[-4:], [
            self.handler._url, '--comment', '-b', self.payload,
        ])

    def test_issue_label_is_literal(self):
        self.handler._pr = False
        self.handler._url = 'https://github.com/example/project/issues/42'
        self.handler.edit(self.payload, '--add-label')
        self.assertEqual(self.arguments()[2:], [
            'issue', 'edit', self.handler._url, '--add-label', self.payload,
        ])

    def test_review_flag_and_ready_have_no_empty_arguments(self):
        self.handler.review('--approve')
        self.assertEqual(self.arguments()[-2:], [self.handler._url, '--approve'])
        self.handler.mark_as_ready()
        self.assertEqual(self.arguments()[-2:], ['ready', self.handler._url])

    def test_rest_fields_are_literal(self):
        endpoint = '/repos/example/project/pulls/42/reviews/7/dismissals'
        self.handler.rest_api(endpoint, method='PUT', query={
            'message': self.payload, 'event': 'DISMISS',
        })
        self.assertEqual(self.arguments(), [
            'api', '-X', 'PUT', '-H', 'Accept: application/vnd.github+json',
            endpoint, '-f', 'message=' + self.payload, '-f', 'event=DISMISS',
        ])

    def test_view_parses_json(self):
        os.environ['FAKE_GH_RESPONSE'] = json.dumps({'labels': [{'name': self.payload}]})
        self.assertEqual(self.handler.view('labels'), [{'name': self.payload}])
        self.assertEqual(self.arguments(), ['pr', 'view', self.handler._url, '--json', 'labels'])

    def test_pagination_preserves_get_method(self):
        os.environ['FAKE_GH_RESPONSE'] = '[{"id": 1}]'
        self.assertEqual(self.handler.query_multi_pages('/repos/example/project/issues/comments'), [{'id': 1}])
        self.assertEqual(self.arguments(), [
            'api', '-X', 'GET', '-H', 'Accept: application/vnd.github+json',
            '/repos/example/project/issues/comments', '-f', 'per_page=100', '-f', 'page=1',
        ])

    def test_command_failure_is_propagated(self):
        os.environ['FAKE_GH_EXIT'] = '1'
        with self.assertRaises(subprocess.CalledProcessError):
            self.handler.mark_as_ready()

    def test_event_arguments_preserve_empty_fields(self):
        args = ['sync_labels.py', 'opened', self.handler._url, 'contributor', '', '', '--info']
        with patch.object(sys, 'argv', args), patch.object(sync_labels, 'GhLabelSynchronizer') as handler:
            sync_labels.main()
        handler.assert_called_once_with(self.handler._url, 'contributor')
        handler.return_value.run.assert_called_once_with(sync_labels.Action.opened, label='', rev_state='')

    def test_scheduled_cleanup(self):
        url = 'https://github.com/example/project'
        with patch.object(sys, 'argv', ['sync_labels.py', url, '--info']), \
                patch.object(sync_labels, 'GhLabelSynchronizer') as handler:
            sync_labels.main()
        handler.assert_called_once_with(url, 'sagetrac-github-bot')
        handler.return_value.run.assert_not_called()


if __name__ == '__main__':
    unittest.main()
