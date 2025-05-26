document.addEventListener('DOMContentLoaded', () => {
  const links = Array.from(
    document.querySelectorAll('.on-this-page a.reference.internal')
  );
  const sections = links
    .map((link) => {
      const href = link.getAttribute('href');
      // Map “#” to top-of-page
      const sec =
        href === '#'
          ? document.documentElement
          : document.getElementById(href.slice(1));
      return sec ? { link, sec } : null;
    })
    .filter((x) => x);

  const lastLink = sections[sections.length - 1].link;

  const onScroll = () => {
    const scrollY = window.scrollY || window.pageYOffset;

    // At very top: only highlight the “#” link
    if (scrollY === 0) {
      links.forEach((l) =>
        l.classList.toggle('active', l.getAttribute('href') === '#')
      );
      return;
    }

    const viewportBottom = scrollY + window.innerHeight;
    const docHeight = document.documentElement.scrollHeight;

    // At bottom: highlight last
    if (viewportBottom >= docHeight - 1) {
      links.forEach((l) => l.classList.toggle('active', l === lastLink));
      return;
    }

    // Normal 15%-down logic
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
  onScroll(); // initialize on load
});
