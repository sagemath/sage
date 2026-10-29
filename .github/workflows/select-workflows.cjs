// Event and path filters for the reusable workflows called by release.yml.
// Run locally with: node --test build/test/ci/select-workflows.test.cjs
const workflows = [
  'build', 'distro', 'meson', 'minimal', 'optional', 'uv', 'doc_html',
  'doc_pdf', 'livedoc', 'homebrew', 'docker', 'release', 'types',
];
const regular = ['build', 'meson', 'minimal', 'homebrew', 'types'];
const docs = ['doc_html', 'doc_pdf'];

module.exports = async ({github, context, core}) => {
  const selected = new Set();
  const add = (...names) => names.forEach(name => selected.add(name));
  const {eventName, payload, ref} = context;
  if (eventName === 'workflow_dispatch') {
    const {target = 'standard', packages = ''} = payload.inputs;
    if (target === 'standard') {
      add(...regular, ...docs, 'distro', 'uv');
    } else if (!workflows.includes(target)) {
      throw new Error(`Unknown CI target: ${target}`);
    } else if (target === 'optional' && !packages.trim()) {
      throw new Error('The optional target requires a nonempty packages input');
    } else {
      add(target);
    }
  } else if (eventName === 'pull_request') {
    add(...regular, ...docs);
    const files = await github.paginate(github.rest.pulls.listFiles, {
      ...context.repo, pull_number: payload.pull_request.number, per_page: 100,
    });
    // Include both sides of renames. At the API's 3000-file limit, run all
    // path-filtered workflows rather than silently miss a relevant change.
    const paths = files.flatMap(f => [f.filename, f.previous_filename].filter(Boolean));
    const truncated = files.length >= 3000;
    if (truncated || paths.some(p =>
      p === '.github/workflows/ci-distro.yml' || p === 'configure.ac' ||
      p.startsWith('build/pkgs/') || p.startsWith('build/make/'))) {
      add('distro');
    }
    if (truncated || paths.some(p =>
      p.endsWith('.build') || p === 'uv.lock' || p.startsWith('subprojects/') ||
      p === '.github/workflows/ci.yml' || /^build\/pkgs\/[^/]+\/distros\/[^/]+\.txt$/.test(p))) {
      add('uv');
    }
    if (truncated || paths.some(p =>
      p === '.github/workflows/push_to_docker_hub.yml' || p === 'docker/Dockerfile')) {
      add('docker');
    }
    if (truncated || paths.includes('.github/workflows/release.yml')) {
      add('release');
    }
  } else if (eventName === 'merge_group') {
    add('build', 'types', ...docs);
  } else if (eventName === 'schedule') {
    add('docker');
  } else if (eventName === 'push') {
    if (ref === 'refs/heads/master' || ref === 'refs/heads/develop') {
      add(...regular);
      if (ref === 'refs/heads/develop') add(...docs);
    } else if (ref.startsWith('refs/tags/')) {
      const tag = ref.slice('refs/tags/'.length);
      // The previous '*' tag filter excludes slashes.
      if (!tag.includes('/')) add('distro', 'uv');
      if (/^\d+\.\d+(?:\.\d+)?(?:\.(?:beta|rc)\d+)?$/.test(tag)) {
        add('livedoc', 'docker', 'release');
      }
    }
  }
  for (const name of workflows) core.setOutput(name, selected.has(name));
};
