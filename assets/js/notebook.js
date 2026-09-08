(() => {
  'use strict';
  const contents = document.querySelector('.notebook-contents');
  if (!contents) return;
  const links = [...contents.querySelectorAll('nav a')];
  const sections = links.map(link => document.getElementById(decodeURIComponent(link.hash.slice(1))));
  const narrow = window.matchMedia('(max-width: 1050px)');
  const setLayout = () => { contents.open = !narrow.matches; };
  setLayout();
  narrow.addEventListener('change', setLayout);
  let active = -1;
  let scheduled = false;
  const update = () => {
    scheduled = false;
    let index = 0;
    sections.forEach((heading, i) => {
      if (heading && heading.getBoundingClientRect().top <= 100) index = i;
    });
    if (index === active) return;
    active = index;
    links.forEach((link, i) => {
      if (i === index) link.setAttribute('aria-current', 'location');
      else link.removeAttribute('aria-current');
    });
  };
  const schedule = () => {
    if (!scheduled) { scheduled = true; window.requestAnimationFrame(update); }
  };
  window.addEventListener('scroll', schedule, { passive: true });
  window.addEventListener('resize', schedule);
  window.addEventListener('load', schedule);
  contents.addEventListener('click', event => {
    if (narrow.matches && event.target.closest('a')) contents.open = false;
  });
  update();
})();
