#!/bin/sh
if [ $# != 2 ]; then
    echo >&2 "usage: $0 BASE_DOC_COMMIT DOC_REPOSITORY"
    echo >&2 "creates CHANGES.html in the docs subdirectory"
    echo >&2 "for the diffs of DOC_REPOSITORY against BASE_DOC_COMMIT"
    exit 1
fi
BASE_DOC_COMMIT="$1"
DOC_REPOSITORY="$2"

mkdir -p ./docs
# Wipe out chronic diffs between old doc and new doc
(cd $DOC_REPOSITORY && find . -name "*.html" | xargs sed -i -e '\;<script type="application/vnd\.jupyter\.widget-state+json">;,\;</script>; d')
# Create CHANGES.html
echo '<html>' > ./docs/CHANGES.html
echo '<head>' >> ./docs/CHANGES.html
echo '<link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/highlight.js/11.9.0/styles/default.min.css">' >> ./docs/CHANGES.html
echo '<script src="https://cdnjs.cloudflare.com/ajax/libs/highlight.js/11.9.0/highlight.min.js"></script>' >> ./docs/CHANGES.html
echo '<script>hljs.highlightAll();</script>' >> ./docs/CHANGES.html
cat >> ./docs/CHANGES.html << EOF
<script>
document.addEventListener('DOMContentLoaded', () => {
const diffSite = 'https://pianomister.github.io/diffsite'
const baseDocURL = 'https://sagemath-tobias.netlify.app'
const diffParagraphs = document.querySelectorAll('p.diff');
diffParagraphs.forEach(paragraph => {
  const rootURL = window.location.origin;
  const docAnchor = paragraph.querySelector('a');  // first "a" element
  const url = new URL(docAnchor.href);
  const path = url.pathname;
  const anchor = document.createElement('a');
  anchor.href = diffSite + '/?url1=' + rootURL + path + '&url2=' + baseDocURL + path;
  anchor.textContent = 'compare with the base';
  anchor.setAttribute('target', '_blank');
  paragraph.appendChild(anchor);
  paragraph.innerHTML += '&nbsp;';
  const hunkAnchors = paragraph.querySelectorAll('a.hunk');
  hunkAnchors.forEach(hunkAnchor => {
    const url = new URL(hunkAnchor.href);
    const path = url.pathname;
    const pathHash = path + url.hash.replace('#', '%23');
    const anchor = document.createElement('a');
    anchor.href = diffSite + '/?url1=' + rootURL + pathHash + '&url2=' + baseDocURL + path;
    anchor.textContent = hunkAnchor.textContent;
    anchor.setAttribute('target', '_blank');
    paragraph.appendChild(anchor);
    paragraph.innerHTML += '&nbsp;';
  });
});
});
</script>
EOF
echo '</head>' >> ./docs/CHANGES.html
echo '<body>' >> ./docs/CHANGES.html
(cd $DOC_REPOSITORY && git diff $BASE_DOC_COMMIT -- *.html; rm -rf .git) > ./docs/diff.txt
/sage/sage -python - << EOF
import os, re, html
with open('./docs/diff.txt', 'r') as f:
    diff_text = f.read()
diff_blocks = re.split(r'^(?=diff --git)', diff_text, flags=re.MULTILINE)
out_blocks = []
for block in diff_blocks:
    match = re.search(r'^diff --git a/(.*) b/\1', block, flags=re.MULTILINE)
    if match:
        path = 'html/' + match.group(1)
        file_path = os.path.join('$DOC_REPOSITORY/..', path)
        with open(file_path, 'r') as file:
            content = file.readlines()
        count = 0
        for line in block.splitlines():
            if line.startswith('@@ -'):
                line_number = int(re.search(r'@@ -(\d+)', line).group(1))
                for i in range(line_number, -1, -1):
                    if content[i].startswith('<'):
                        count += 1
                        content[i] = f'<span id="hunk{count}" style="visibility: hidden;"></span>' + content[i]
                        break
        with open(file_path, 'w') as file:
            file.writelines(content)
        hunks = '&nbsp;'.join(f'<a href="{path}#hunk{i+1}" class="hunk" target="_blank">#{i+1}</a>' for i in range(count))
        out_blocks.append(f'<p class="diff"><a href="{path}">{path}</a>&nbsp;' + hunks + '&emsp;</p>'
                            + '\n<pre><code class="language-diff">'
                            + html.escape(block).strip() + '</code></pre>')
output_text = '\n'.join(out_blocks)
with open('./docs/diff.html', 'w') as f:
    f.write(output_text)
EOF
cat ./docs/diff.html >> ./docs/CHANGES.html
echo '</body>' >> ./docs/CHANGES.html
echo '</html>' >>./docs/CHANGES.html
rm ./docs/diff.txt ./docs/diff.html
