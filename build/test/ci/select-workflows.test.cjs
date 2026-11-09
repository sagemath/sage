const assert = require('node:assert/strict');
const {test} = require('node:test');
const select = require('../../../.github/workflows/select-workflows.cjs');
const regular = ['build', 'meson', 'minimal', 'homebrew', 'types'];
const docs = ['doc_html', 'doc_pdf'];

async function selection(eventName, {paths = [], ref = '', inputs = {}} = {}) {
  const output = {};
  await select({
    context: {
      eventName, ref, repo: {owner: 'sagemath', repo: 'sage'},
      payload: {inputs, pull_request: {number: 1}},
    },
    github: {
      rest: {pulls: {listFiles: 'listFiles'}},
      paginate: async (method, args) => {
        assert.equal(method, 'listFiles');
        assert.equal(args.pull_number, 1);
        return paths.map(p => typeof p === 'string' ? {filename: p} : p);
      },
    },
    core: {setOutput: (key, value) => { output[key] = value; }},
  });
  assert.equal(Object.keys(output).length, 13);
  return Object.keys(output).filter(key => output[key]).sort();
}

const cases = [
  ['ordinary PR', 'pull_request', {}, [...regular, ...docs]],
  ['distro change', 'pull_request', {paths: ['build/pkgs/flint/checksums.ini']}, [...regular, ...docs, 'distro']],
  ['Meson change', 'pull_request', {paths: ['src/sage/meson.build']}, [...regular, ...docs, 'uv']],
  ['distro mapping', 'pull_request', {paths: ['build/pkgs/flint/distros/debian.txt']}, [...regular, ...docs, 'distro', 'uv']],
  ['release workflow', 'pull_request', {paths: ['.github/workflows/release.yml']}, [...regular, ...docs, 'release']],
  ['Dockerfile', 'pull_request', {paths: ['docker/Dockerfile']}, [...regular, ...docs, 'docker']],
  ['renamed file', 'pull_request', {paths: [{filename: 'old-configure.ac', previous_filename: 'configure.ac'}]}, [...regular, ...docs, 'distro']],
  ['merge queue', 'merge_group', {}, ['build', 'types', ...docs]],
  ['develop', 'push', {ref: 'refs/heads/develop'}, [...regular, ...docs]],
  ['master', 'push', {ref: 'refs/heads/master'}, regular],
  ['other branch', 'push', {ref: 'refs/heads/topic'}, []],
  ['ordinary tag', 'push', {ref: 'refs/tags/test'}, ['distro', 'uv']],
  ['tag with slash', 'push', {ref: 'refs/tags/test/1'}, []],
  ['scheduled Docker build', 'schedule', {}, ['docker']],
  ['manual standard', 'workflow_dispatch', {}, [...regular, ...docs, 'distro', 'uv']],
];
for (const tag of ['11.0', '11.0.1', '11.0.beta5', '11.0.rc1', '11.0.1.beta2', '11.0.1.rc2']) {
  cases.push([`release ${tag}`, 'push', {ref: `refs/tags/${tag}`}, ['distro', 'uv', 'livedoc', 'docker', 'release']]);
}
for (const target of [...regular, ...docs, 'distro', 'uv', 'livedoc', 'docker', 'release', 'optional']) {
  cases.push([`manual ${target}`, 'workflow_dispatch', {inputs: {target, packages: 'foo,bar'}}, [target]]);
}
for (const [name, event, options, expected] of cases) {
  test(name, async () => assert.deepEqual(await selection(event, options), [...expected].sort()));
}
test('reject missing optional packages', async () => {
  await assert.rejects(selection('workflow_dispatch', {inputs: {target: 'optional'}}), /requires a nonempty packages/);
});
test('reject unknown manual target', async () => {
  await assert.rejects(selection('workflow_dispatch', {inputs: {target: 'typo'}}), /Unknown CI target/);
});
test('conservative selection at the API limit', async () => {
  assert.deepEqual(await selection('pull_request', {paths: Array(3000).fill('README.md')}),
    [...regular, ...docs, 'distro', 'uv', 'docker', 'release'].sort());
});
