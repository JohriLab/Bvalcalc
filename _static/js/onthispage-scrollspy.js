document.addEventListener('DOMContentLoaded', () => {
  const links = Array.from(
    document.querySelectorAll('.on-this-page a.reference.internal')
  );
  const sections = links
    .map((link) => {
      const href = link.getAttribute('href');
      const sec =
        href === '#'
          ? document.documentElement
          : document.getElementById(href.slice(1));
      return sec ? { link, sec } : null;
    })
    .filter((x) => x);

  const lastLink = sections[sections.length - 1].link;

  let forcedLink = null;
  let manualScroll = false;

  // Clicking a TOC link forces that highlight until the first manual scroll
  links.forEach((link) => {
    link.addEventListener('click', () => {
      forcedLink = link;
      manualScroll = false; // reset our “I’ve scrolled” flag
      links.forEach((l) => l.classList.toggle('active', l === forcedLink));
    });
  });

  // Flag real user scrolls (wheel or touch)
  window.addEventListener(
    'wheel',
    () => {
      manualScroll = true;
    },
    { passive: true }
  );
  window.addEventListener(
    'touchmove',
    () => {
      manualScroll = true;
    },
    { passive: true }
  );

  const onScroll = () => {
    // If we’ve forced one but haven’t yet manually scrolled, bail out
    if (forcedLink && !manualScroll) {
      return;
    }

    // Once you *do* manually scroll, we clear the forced state
    if (forcedLink && manualScroll) {
      forcedLink = null;
    }

    const scrollY = window.scrollY || window.pageYOffset;

    // 1) Top of page → only “#” link
    if (scrollY === 0) {
      links.forEach((l) =>
        l.classList.toggle('active', l.getAttribute('href') === '#')
      );
      return;
    }

    const viewportBottom = scrollY + window.innerHeight;
    const docHeight = document.documentElement.scrollHeight;

    // 2) Bottom of page → last link
    if (viewportBottom >= docHeight - 1) {
      links.forEach((l) => l.classList.toggle('active', l === lastLink));
      return;
    }

    // 3) Otherwise → your 15%–down trigger
    const triggerPoint = scrollY + window.innerHeight * 0.15;
    let current = null;
    for (const { link, sec } of sections) {
      if (sec.offsetTop <= triggerPoint) {
        current = link;
      }
    }

    if (current) {
      links.forEach((l) => l.classList.toggle('active', l === current));
    } else {
      links.forEach((l) => l.classList.remove('active'));
    }
  };

  window.addEventListener('scroll', onScroll);
  window.addEventListener('hashchange', onScroll);
  onScroll(); // initial run

  // ——— Scroll-to-bottom button handler ———
  const scrollBtn = document.querySelector('.scroll-to-bottom-btn');
  const scroller = document.querySelector('.wy-side-scroll');
  if (scrollBtn && scroller) {
    scrollBtn.addEventListener('click', (e) => {
      e.preventDefault();
      scroller.scrollTo({
        top: scroller.scrollHeight,
        behavior: 'smooth',
      });
    });
  }
});
