// Exercise the actual navigation script without requiring a browser.
const fs = require('node:fs');
const vm = require('node:vm');
const assert = require('node:assert/strict');
const path = require('node:path');
const script = fs.readFileSync(path.join(__dirname, '../assets/js/notebook.js'), 'utf8');
const ids = ['Overview', 'Least squares', 'Écart 10%'];
let positions = [0, 400, 900];
const links = ids.map(id => ({hash: '#' + encodeURIComponent(id), attrs: {},
  setAttribute(key, value) {this.attrs[key] = value;},
  removeAttribute(key) {delete this.attrs[key];}}));
const sections = Object.fromEntries(ids.map((id, i) => [id, {getBoundingClientRect: () => ({top: positions[i]})}]));
const events = {};
const contents = {open: true, querySelectorAll: () => links,
  addEventListener: (name, handler) => {events[name] = handler;}};
const media = {matches: false, addEventListener: (_, handler) => {media.change = handler;}};
const frame = [];
vm.runInNewContext(script, {document: {querySelector: () => contents, getElementById: id => sections[id]},
  window: {matchMedia: () => media, addEventListener: (name, handler) => {events[name] = handler;},
    requestAnimationFrame: handler => frame.push(handler)}});
assert.equal(contents.open, true);
assert.equal(links[0].attrs['aria-current'], 'location');
positions = [-500, -100, 400]; events.scroll(); events.scroll();
assert.equal(frame.length, 1, 'Scroll updates should be batched');
frame.shift()();
assert.equal(links[1].attrs['aria-current'], 'location');
assert.equal(links[0].attrs['aria-current'], undefined);
positions = [-1100, -700, 20]; events.scroll(); frame.shift()();
assert.equal(links[2].attrs['aria-current'], 'location');
media.matches = true; media.change(); assert.equal(contents.open, false);
contents.open = true; events.click({target: {closest: () => links[1]}});
assert.equal(contents.open, false, 'Mobile navigation closes after choosing a section');
media.matches = false; media.change(); assert.equal(contents.open, true);
console.log('Notebook navigation passed: current section, scroll batching, Unicode anchors, responsive menu, and mobile selection.');
