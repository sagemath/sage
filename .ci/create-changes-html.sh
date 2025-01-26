#!/bin/sh
if [ $# != 2 ]; then
    echo >&2 "Usage: $0 DIFF_TEXT DOC_REPO"
    echo >&2 "This script generates a CHANGES.html file in the current directory"
    echo >&2 "and adds anchor targets in the documents within DOC_REPO"
    echo >&2 "based on the diff hunks in the DIFF_TEXT file."
    exit 1
fi
DIFF_TEXT="$1"
DOC_REPOSITORY="$2"

# Create CHANGES.html
echo '<html>' > CHANGES.html
echo '<head>' >> CHANGES.html
echo '<link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/highlight.js/11.9.0/styles/default.min.css">' >> CHANGES.html
echo '<script src="https://cdnjs.cloudflare.com/ajax/libs/highlight.js/11.9.0/highlight.min.js"></script>' >> CHANGES.html
echo '<script>hljs.highlightAll();</script>' >> CHANGES.html
cat >> CHANGES.html << EOF
<style>
  p.diff a:first-child {
    font-weight: bold;
    font-size: x-large;
  }
</style>
<script>
document.addEventListener('DOMContentLoaded', () => {
// This URL is hardcoded in the file .github/workflows/doc-publish.yml.
// See NETLIFY_ALIAS of the "Deploy to Netlify" step.
const baseDocURL = 'https://doc-develop--sagemath.netlify.app'
const diffSite = 'https://pianomister.github.io/diffsite'
const diffParagraphs = document.querySelectorAll('p.diff');
diffParagraphs.forEach(paragraph => {
  const rootURL = window.location.origin;
  const docAnchor = paragraph.querySelector('a');
  const url = new URL(docAnchor.href);
  const path = url.pathname;
  const anchor = document.createElement('a');
  anchor.href = diffSite + '/?url1=' + rootURL + path + '&url2=' + baseDocURL + path;
  anchor.textContent = 'compare with the base';
  anchor.setAttribute('target', '_blank');
  paragraph.innerHTML += '&nbsp;&nbsp;';
  paragraph.appendChild(anchor);
  const hunks = paragraph.parentNode.querySelectorAll('p.hunk');
  hunks.forEach(hunk => {
    const hunkAnchor = hunk.querySelector('a');
    const url = new URL(hunkAnchor.href);
    const path = url.pathname;
    const pathHash = path + url.hash.replace('#', '%23');
    const anchor = document.createElement('a');
    anchor.href = diffSite + '/?url1=' + rootURL + pathHash + '&url2=' + baseDocURL + path;
    anchor.textContent = 'compare with the base';
    anchor.setAttribute('target', '_blank');
    hunk.innerHTML += '&nbsp;&nbsp;';
    hunk.appendChild(anchor);
  });
});
});
</script>
EOF
echo '</head>' >> CHANGES.html
echo '<body>' >> CHANGES.html
python3 - << EOF
import os, re, html
from itertools import chain
with open('$DIFF_TEXT', 'r') as f:
    diff_text = f.read()
diff_blocks = re.split(r'^(?=diff --git)', diff_text, flags=re.MULTILINE)
out_blocks = []
for block in diff_blocks:
    match = re.search(r'^diff --git a/(.*) b/\1', block, flags=re.MULTILINE)
    if match:
        path = match.group(1)
        file_path = os.path.join('$DOC_REPOSITORY', path)
        try:
            with open(file_path, 'r') as file:
                content = file.readlines()
        except FileNotFoundError:
            content = []
        count = 0
        hunks = []
        hunk_lines = []
        in_hunk = False
        for line in block.splitlines():
            if line.startswith('@@ -'):
                if hunk_lines:
                    hunks.append('<pre><code class="language-diff">'
                                 + html.escape('\n'.join(hunk_lines)).strip()
                                 + '</code></pre>')
                    hunk_lines = []
                search_result = re.search(r'@@ -(\d+),(\d+) \+(\d+),(\d+)', line)
                if search_result:
                    line_number = int(search_result.group(3)) - 1
                    span = int(search_result.group(4))
                    for i in chain(range(line_number, line_number + span), range(line_number - 1, -1, -1)):
                        try:
                            ln = content[i]
                        except IndexError:
                            continue
                        for idx, char in enumerate(ln):
                            if not char.isspace():
                                break
                        else:
                            idx = len(ln)
                        if ln.startswith('<', idx) and not ln.startswith('</', idx):
                            count += 1
                            content[i] = ln[:idx] + f'<span id="hunk{count}" style="visibility: hidden;"></span>' + ln[idx:]
                            hunks.append(f'<p class="hunk"><a href="{path}#hunk{count}" class="hunk" target="_blank">hunk #{count}</a></p>')
                            break
            hunk_lines.append(line)
        if hunk_lines:
            hunks.append('<pre><code class="language-diff">'
                          + html.escape('\n'.join(hunk_lines)).strip()
                          + '</code></pre>')
        if content:
            with open(file_path, 'w') as file:
                file.writelines(content)
        out_blocks.append(f'<div class="diff"><p class="diff"><a href="{path}">{path}</a></p>\n' + '\n'.join(hunks) + '\n</div>')
output_text = '\n'.join(out_blocks)
with open('diff.html', 'w') as f:
    f.write(output_text)
EOF
cat diff.html >> CHANGES.html
echo '</body>' >> CHANGES.html
echo '</html>' >> CHANGES.html
rm diff.html
