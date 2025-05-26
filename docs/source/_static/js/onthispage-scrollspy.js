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

  const onScroll = () => {
    const scrollY = window.scrollY || window.pageYOffset;
    // Calculate the 70% down-the-screen trigger point
    const triggerPoint = scrollY + window.innerHeight * 0.15;
    let current = sections[0].link;

    for (const { link, sec } of sections) {
      if (sec.offsetTop <= triggerPoint) {
        current = link;
      }
    }

    links.forEach((l) => l.classList.toggle('active', l === current));
  };

  window.addEventListener('scroll', onScroll);
  onScroll(); // initialize on load
});
