/* Search is entirely local: no queries or browsing data leave the page. */
(() => {
  const root = document.querySelector('[data-catalog]');
  if (!root) return;
  const input = root.querySelector('#tour-search');
  const buttons = [...root.querySelectorAll('[data-category]')].filter(el => el.tagName === 'BUTTON');
  const cards = [...root.querySelectorAll('[data-tour]')];
  const count = root.querySelector('#result-count');
  const empty = root.querySelector('#empty-results');
  const normalize = value => value.normalize('NFD').replace(/[\u0300-\u036f]/g, '').toLowerCase().replace(/ℓ/g, 'l').replace(/¹/g, '1');
  const index = new Map(cards.map(card => [card, normalize(card.dataset.search)]));
  let category = 'all';
  function render(updateUrl = true) {
    const terms = normalize(input.value.trim()).split(/\s+/).filter(Boolean);
    let visible = 0;
    cards.forEach(card => {
      const matches = (category === 'all' || card.dataset.category === category) && terms.every(term => index.get(card).includes(term));
      card.hidden = !matches;
      visible += Number(matches);
    });
    buttons.forEach(button => button.setAttribute('aria-pressed', String(button.dataset.category === category)));
    count.textContent = `${visible} ${visible === 1 ? 'tour' : 'tours'}${visible === cards.length ? '' : ` of ${cards.length}`}`;
    empty.hidden = visible !== 0;
    if (updateUrl) {
      const url = new URL(location.href);
      input.value.trim() ? url.searchParams.set('q', input.value.trim()) : url.searchParams.delete('q');
      category === 'all' ? url.searchParams.delete('topic') : url.searchParams.set('topic', category);
      history.replaceState(null, '', url);
    }
  }
  function fromUrl() {
    const params = new URLSearchParams(location.search);
    input.value = params.get('q') || '';
    category = buttons.some(button => button.dataset.category === params.get('topic')) ? params.get('topic') : 'all';
    render(false);
  }
  input.addEventListener('input', () => render());
  root.querySelector('form').addEventListener('submit', event => { event.preventDefault(); render(); });
  buttons.forEach(button => button.addEventListener('click', () => { category = button.dataset.category; render(); }));
  root.querySelector('#reset-search').addEventListener('click', () => { input.value = ''; category = 'all'; render(); input.focus(); });
  window.addEventListener('popstate', fromUrl);
  fromUrl();
})();
