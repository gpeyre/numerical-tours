// Run the real search handlers against a small DOM fixture, without a browser.
const fs = require('node:fs');
const vm = require('node:vm');
const assert = require('node:assert/strict');
const path = require('node:path');
const site = path.resolve(__dirname, '..');
const data = JSON.parse(fs.readFileSync(path.join(site, '_data/python_catalog.json'), 'utf8'));
function element(dataset = {}, tagName = 'DIV') {
  return {dataset, tagName, hidden: false, attrs: {}, events: {}, value: '',
    setAttribute(k, v) {this.attrs[k] = v;},
    addEventListener(k, fn) {this.events[k] = fn;},
    focus() {this.focused = true;}};
}
const cards = data.map(t => element({search: t.search, category: t.category}));
const buttons = ['all', ...new Set(data.map(t => t.category))].map(category => element({category}, 'BUTTON'));
const input = element(), count = element(), empty = element(), form = element(), reset = element();
const elements = {'#tour-search': input, '#result-count': count, '#empty-results': empty, form, '#reset-search': reset};
const root = {querySelector: s => elements[s], querySelectorAll: s => s === '[data-tour]' ? cards : buttons};
const location = {href: 'https://example.com/python/?q=Sinkhorn', search: '?q=Sinkhorn'};
const history = {replaceState(_, __, url) {location.href = String(url); location.search = url.search;}};
const events = {};
const context = {document: {querySelector: () => root}, location, history, URL, URLSearchParams,
  window: {addEventListener(k, fn) {events[k] = fn;}}};
vm.runInNewContext(fs.readFileSync(path.join(site, 'assets/js/catalog.js'), 'utf8'), context);
const visible = () => data.filter((t, i) => !cards[i].hidden);
assert(visible().length > 0);
assert(visible().every(t => t.search.toLowerCase().includes('sinkhorn')));
input.value = 'francais'; input.events.input();
assert(visible().some(t => t.slug === 'introduction_6_elementary_fr'));
input.value = 'neural network'; input.events.input();
assert(visible().length >= 4);
buttons.find(b => b.dataset.category === 'Optimal transport').events.click();
assert.equal(visible().length, 0);
assert.equal(empty.hidden, false);
reset.events.click();
assert.equal(visible().length, 58);
assert.equal(input.focused, true);
assert.equal(location.search, '');
input.value = '<script>alert(1)</script>'; input.events.input();
assert.equal(visible().length, 0);
location.search = '?topic=not-a-category&q=wavelet'; events.popstate();
assert(visible().length > 0);
assert.equal(buttons[0].attrs['aria-pressed'], 'true');
assert(count.textContent.includes('of 58'));
form.events.submit({preventDefault() {}});
console.log('Search checks passed: query, accents, multiple words, filters, empty state, reset, URL state, and literal input.');
