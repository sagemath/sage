// Change the editor theme of all CodeMirror cells according to the furo (dark) mode
function changeTheme(editor, theme) {
  if (theme === 'dark') {
    editor.setOption('theme', 'monokai'); // the same with pygments dark style in conf.py
  }
  else {
    editor.setOption('theme', 'default');
  }
}

// Uses the theme data of the document.body element set by setTheme function
// defined in https://github.com/pradyunsg/furo/blob/main/src/furo/assets/scripts/furo.js
const body = document.body;
const observer1 = new MutationObserver((mutationsList) => {
  for (let mutation of mutationsList) {
    if (mutation.type === 'attributes' && mutation.attributeName === 'data-theme') {
      const theme = body.dataset.theme;
      const querySet = document.querySelectorAll('.CodeMirror');
      for (var i = 0; i < querySet.length; i++) {
        changeTheme(querySet[i].CodeMirror, theme);
      }}}});
observer1.observe(body, { attributes: true });

// Change the editor theme of a new CodeMirror cell.
const callback = function(mutationsList, observer) {
  for(const mutation of mutationsList) {
    if (mutation.type === 'childList') {
      const theme = body.dataset.theme;
      for (const addedNode of mutation.addedNodes) {
        if (addedNode.classList && addedNode.classList.contains('CodeMirror')) {
          changeTheme(addedNode.CodeMirror, theme);
        }}}}};
const observer2 = new MutationObserver(callback);
observer2.observe(document.getElementsByClassName("content")[0], { childList: true, subtree: true });
