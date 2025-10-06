// Change the editor theme according to the furo light/dark/auto mode
function changeTheme(editor, theme) {
  if (theme === 'dark') {
    editor.setOption('theme', 'monokai'); // the same with pygments dark style in conf.py
  } else if (theme === 'auto' && window.matchMedia('(prefers-color-scheme: dark)').matches) {
    editor.setOption('theme', 'monokai');
  } else {
    editor.setOption('theme', 'default');
  }
}

// Change the editor theme of all CodeMirror cells
function changeThemeAll(theme) {
  const querySet = document.querySelectorAll('.CodeMirror');
  for (var i = 0; i < querySet.length; i++) {
    changeTheme(querySet[i].CodeMirror, theme);
  }
}

// Use the theme data of the body element set by setTheme function
// defined in https://github.com/pradyunsg/furo/blob/main/src/furo/assets/scripts/furo.js
const body = document.body;
const observer1 = new MutationObserver((mutationsList) => {
  for (let mutation of mutationsList) {
    if (mutation.type === 'attributes' && mutation.attributeName === 'data-theme') {
      const theme = body.dataset.theme;
      changeThemeAll(theme);
    }
  }
});
observer1.observe(body, { attributes: true });


// In the furo auto mode, we watch prefers-color-scheme and use the theme data
// of the body element to change the CodeMirror editor theme
const prefersDarkMode = window.matchMedia('(prefers-color-scheme: dark)');

function handlePrefersColorSchemeChange(e) {
  const theme = body.dataset.theme;
  if (theme === 'auto') {
    changeThemeAll(theme);
  }
}

prefersDarkMode.addEventListener('change', handlePrefersColorSchemeChange);


// Change the editor theme of a new CodeMirror cell
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


//
// Version selector
//

function fetchVersions() {
    try {
        let menu = document.getElementById('versions-menu');

        // For the origin of the this site, see .github/workflows/doc-publish.yml
        fetch('https://doc-release--sagemath.netlify.app/html/en/versions.txt')
            .then(response => {
                if (!response.ok) {
                    throw new Error('Network response was not ok ' + response.statusText);
                }
                return response.text();
            })
            .then(text => {
                const lines = text.split('\n');
                lines.forEach(line => {
                    if (!line.startsWith('#')) { // Ignore the comment line
                        let [ver, url] = line.split(' ');
                        if (ver && url) {
                            if (!url.startsWith('https://')) {
                                url = 'https://' + url;
                            }
                            let option = document.createElement('option');
                            option.value = url;
                            option.text = ver;
                            menu.add(option);
                        }
                    }
                });
            })
            .catch(error => {
                console.log('Failed to fetch versions.txt file.');
            });
    } catch (error) {
        console.log('Failed to fetch versions.txt file.');
    }
}

fetchVersions()

// Function to change the version based on versions menu selection
function changeVersion() {
    let select_element = document.getElementById("versions-menu");
    let selected_ver = select_element.options[select_element.selectedIndex].text;
    let selected_url = select_element.value;
    if (selected_url) {
        if (window.location.protocol == 'file:') {
            let pathname = window.location.pathname;
            let cutoff_point = pathname.indexOf('/html');
            if (cutoff_point !== -1) {
                pathname = pathname.substring(cutoff_point);
                window.location.href = selected_url + pathname;
            } else {
                window.location.href = selected_url + '/index.html';
            }
        } else {
            window.location.href = selected_url + window.location.pathname;
        }
    }
}


// Listen to the kernel status changes
// https://thebe.readthedocs.io/en/stable/events.html
thebelab.on("status", function (evt, data) {
  if (data.status === 'building') {
    const elements = document.querySelectorAll('.thebelab-cell');
    elements.forEach(element => {
      element.style.filter = 'opacity(50%)';
    });
    const element = document.getElementById("thebelab-activate-button");
    element.innerHTML = "Building";
    element.style.right = '.4rem';
  }
  else if (data.status === 'built') {
    const elements = document.querySelectorAll('.thebelab-cell');
    elements.forEach(element => {
      element.style.filter = 'opacity(60%)';
    });
    const element = document.getElementById("thebelab-activate-button");
    element.innerHTML = "Built";
    element.style.right = '.4rem';
  }
  else if (data.status === 'launching') {
    const elements = document.querySelectorAll('.thebelab-cell');
    elements.forEach(element => {
      element.style.filter = 'opacity(70%)';
    });
    const element = document.getElementById("thebelab-activate-button");
    element.innerHTML = "Launching";
    element.style.right = '.4rem';
  }
  else if (data.status === 'failed') {
    const elements = document.querySelectorAll('.thebelab-cell');
    elements.forEach(element => {
      element.style.filter = 'opacity(50%)';
    });
    const element = document.getElementById("thebelab-activate-button");
    element.innerHTML = 'Failed: ' + data.message;
    element.style.right = '.4rem';
    element.style.width = 'auto';
    element.style.color = 'red';
  }
  else if (data.status === 'ready') {
    const elements = document.querySelectorAll('.thebelab-cell');
    elements.forEach(element => {
      element.style.filter = 'opacity(100%)';
    });
    const element = document.getElementById("thebelab-activate-button");
    element.innerHTML = "Ready";
    element.style.right = null;
    // Run custom code when the kernel is ready
    const kernel = data.kernel;
    kernel.requestExecute({code: "%display latex"});
  }
});


// Activate Thebe when "Sage Live" tab is clicked
document.querySelectorAll('input[class="tab-input"]').forEach((elem) => {
    elem.addEventListener("click", function(event) {
        if (elem.nextElementSibling) {
            if (elem.nextElementSibling.nextElementSibling) {
                if (elem.nextElementSibling.nextElementSibling.querySelector('div[class="thebelab-code"]')) {
                    initThebelab();
                }
            }
        }
    });
});
