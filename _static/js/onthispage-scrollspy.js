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

  // Scroll-to-bottom click handler
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

  // Hide button when scroller is at (or within 10px of) its bottom
  if (scroller && scrollBtn) {
    const updateBtn = () => {
      const distanceFromBottom =
        scroller.scrollHeight - scroller.clientHeight - scroller.scrollTop;
      scrollBtn.style.display = distanceFromBottom <= 36 ? 'none' : '';
    };
    updateBtn();
    scroller.addEventListener('scroll', updateBtn, { passive: true });
    window.addEventListener('resize', updateBtn, { passive: true });
  }

  // ─────────────── Extra Footer Injection ───────────────
  const contentInfo = document.querySelector('div[role="contentinfo"]');
  if (contentInfo) {
    // Build the extra text container
    const extra = document.createElement('div');
    extra.className = 'extra-footer-text';
    extra.innerHTML = `
      <p>
        Citation: Marsh JI, Johri, PJ. The B-value calculator: efficient analytical calculation of expected
        genetic diversity with background selection. bioRxiv. doi:
          <a href="https://biorxiv" target="_blank" rel="noopener noreferrer">10.1101/2063.10.10.561727</a>
      </p>
      <div>      
        <div>
        Made with <a href="https://www.sphinx-doc.org/" target="_blank" rel="noopener noreferrer">Sphinx</a>
          using an adapted <a href="https://github.com/rtfd/sphinx_rtd_theme" target="_blank" rel="noopener noreferrer">theme</a> 
          from <a href="https://readthedocs.org/" target="_blank" rel="noopener noreferrer">Read the Docs</a> <br/>
        </div>
        <a href="https://test.pypi.org/project/bvalcalc/" target="_blank" rel="noopener noreferrer">PyPI</a> | 
        <a href="https://github.com/JohriLab/BvalueCalculator/issues" target="_blank" rel="noopener noreferrer">Report a bug</a> | 
        <a href="https://github.com/JohriLab/BvalueCalculator" target="_blank" rel="noopener noreferrer">Source code</a> 
      </div>
    `;

    // Show/remove based on viewport width
    function updateFooterExtras() {
      if (window.innerWidth <= 1304) {
        if (!contentInfo.contains(extra)) contentInfo.appendChild(extra);
      } else {
        if (contentInfo.contains(extra)) contentInfo.removeChild(extra);
      }
    }

    updateFooterExtras();
    window.addEventListener('resize', updateFooterExtras);
  }

  // ── Make the bottom two captions clickable ──
  document
    .querySelectorAll('.wy-menu-vertical > p.caption:nth-last-of-type(-n+2)')
    .forEach((cap) => {
      const nextUL = cap.nextElementSibling;
      if (!nextUL || nextUL.tagName.toLowerCase() !== 'ul') return;

      const childLink = nextUL.querySelector('li > a.reference.internal');
      if (!childLink) return;

      const href = childLink.getAttribute('href');
      const span = cap.querySelector('span.caption-text');
      if (!span) return;

      const text = span.textContent.trim();
      span.textContent = '';

      const a = document.createElement('a');
      a.href = href;
      a.textContent = text;

      // ── copy “current” if the hidden link already has it ──
      if (childLink.classList.contains('current')) {
        a.classList.add('current');
      }

      span.appendChild(a);
    });
  document
    .querySelectorAll('.wy-menu-vertical > p.caption .caption-text a')
    .forEach((a) => {
      const href = a.getAttribute('href');
      if (
        window.location.pathname.endsWith(href) ||
        window.location.href.endsWith(href)
      ) {
        // Force white background even if CSS has higher specificity
        a.style.setProperty('background-color', '#fff', 'important');
      }
    });
});
