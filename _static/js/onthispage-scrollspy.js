document.addEventListener('DOMContentLoaded', () => {
  const links = Array.from(
    document.querySelectorAll('.on-this-page a.reference.internal')
  );
  const sections = links
    .map((link) => {
      const id = link.getAttribute('href').slice(1);
      const sec = document.getElementById(id);
      return sec ? { link, sec } : null;
    })
    .filter((x) => x);

  const lastLink = sections[sections.length - 1].link;

  const onScroll = () => {
    const scrollY = window.scrollY || window.pageYOffset;

    // 1) At the very top, clear all
    if (scrollY === 0) {
      links.forEach((l) => l.classList.remove('active'));
      return;
    }

    const viewportBottom = scrollY + window.innerHeight;
    const docHeight = document.documentElement.scrollHeight;

    // 2) At the very bottom, force-last
    if (viewportBottom >= docHeight - 1) {
      links.forEach((l) => l.classList.toggle('active', l === lastLink));
      return;
    }

    // 3) Only highlight once first section has crossed 15%
    const triggerPoint = scrollY + window.innerHeight * 0.15;
    let current = null;

    for (const { link, sec } of sections) {
      if (sec.offsetTop <= triggerPoint) {
        current = link;
      }
    }

    // Apply or clear
    if (current) {
      links.forEach((l) => l.classList.toggle('active', l === current));
    } else {
      links.forEach((l) => l.classList.remove('active'));
    }
  };

  window.addEventListener('scroll', onScroll);
  onScroll(); // initialize on load
});
