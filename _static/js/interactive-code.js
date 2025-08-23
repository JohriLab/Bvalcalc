// Function to initialize interactive code block
function initializeCodeBlock(codeBlock) {
  codeBlock.contentEditable = true;
  codeBlock.spellcheck = false;
  codeBlock.style.outline = 'none';
  codeBlock.style.cursor = 'text';
  codeBlock.style.caretColor = 'black';
  codeBlock.style.userSelect = 'text';
  codeBlock.style.whiteSpace = 'pre';
  codeBlock.style.paddingTop = '25px';

  // Sanitize input to plain text only
  codeBlock.addEventListener('paste', (e) => {
    e.preventDefault();
    const text = e.clipboardData.getData('text/plain');
    document.execCommand('insertText', false, text);
  });

  // Note: HTML sanitization is handled by the paste event only
  // This prevents cursor jumping issues during normal typing

  // Create copy button
  const copyButton = document.createElement('button');
  copyButton.innerHTML = `
    <svg width="14" height="14" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2">
      <rect x="9" y="9" width="13" height="13" rx="2" ry="2"></rect>
      <path d="M5 15H4a2 2 0 0 1-2-2V4a2 2 0 0 1 2-2h9a2 2 0 0 1 2 2v1"></path>
    </svg>
  `;
  copyButton.style.cssText = `
    position: absolute;
    top: 3px;
    right: 31px;
    padding: 2px;
    background: transparent;
    border: none;
    cursor: pointer;
    display: block;
    z-index: 20;
    color: #FF8C00;
  `;

  // Insert button into the parent container
  const parent = codeBlock.parentElement;
  parent.style.position = 'relative';
  parent.appendChild(copyButton);

  // Create reset button
  const resetButton = document.createElement('button');
  resetButton.innerHTML = `
    <svg width="14" height="14" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2">
      <path d="M3 12a9 9 0 0 1 9-9 9.75 9.75 0 0 1 6.74 2.74L21 8"></path>
      <path d="M21 3v5h-5"></path>
      <path d="M21 12a9 9 0 0 1-9 9 9.75 9.75 0 0 1-6.74-2.74L3 16"></path>
      <path d="M3 21v-5h5"></path>
    </svg>
  `;
  resetButton.style.cssText = `
    position: absolute;
    top: 3px;
    right: 6px;
    padding: 2px;
    background: transparent;
    border: none;
    cursor: pointer;
    display: block;
    z-index: 20;
    color: #FF8C00;
  `;
  parent.appendChild(resetButton);

  // Determine language and add appropriate tag
  const languageTag = document.createElement('div');

  // Check if this is a Python block by looking at the parent's class
  const isPython = parent.parentElement.classList.contains('highlight-python');
  const isBash = parent.parentElement.classList.contains('highlight-bash');

  if (isPython) {
    languageTag.textContent = 'python';
  } else if (isBash) {
    languageTag.textContent = 'bash';
  } else {
    languageTag.textContent = 'code';
  }

  languageTag.style.cssText = `
    position: absolute;
    top: -1px;
    left: 2px;
    background: transparent;
    color: #FF8C00;
    padding: 5px 9px;
    font-size: 11px;
    font-weight: 500;
    z-index: 5;
  `;
  parent.appendChild(languageTag);

  codeBlock.addEventListener('focus', () => {
    codeBlock.style.border = '1px solid black';
  });

  codeBlock.addEventListener('blur', () => {
    codeBlock.style.border = 'none';
  });

  // Create message element
  const message = document.createElement('span');
  message.textContent = 'copied to clipboard!';
  message.style.cssText = `
    position: absolute;
    top: 3px;
    right: 47px;
    color: #495057;
    font-size: 11px;
    opacity: 0;
    pointer-events: none;
    z-index: 15;
    background: transparent;
    padding: 2px 4px;
  `;
  parent.appendChild(message);

  // Store original element
  const originalElement = codeBlock.cloneNode(true);

  // Copy functionality
  copyButton.addEventListener('click', (e) => {
    e.stopPropagation();

    // Get the content and filter out comment lines
    const rawText = codeBlock.innerText || codeBlock.textContent;
    const lines = rawText.split(/\r?\n/);
    const filteredLines = lines.filter((line) => !line.trim().startsWith('#'));
    const cleanText = filteredLines.join('\n').trim();

    // Copy the filtered text to clipboard
    navigator.clipboard.writeText(cleanText);

    // Show message and keep it visible
    message.style.opacity = '1';
  });

  // Reset functionality
  resetButton.addEventListener('click', (e) => {
    e.stopPropagation();
    const newElement = originalElement.cloneNode(true);
    codeBlock.parentNode.replaceChild(newElement, codeBlock);
    message.style.opacity = '0';

    // Reinitialize the new element
    initializeCodeBlock(newElement);
  });
}

document.addEventListener('DOMContentLoaded', () => {
  const codeBlocks = document.querySelectorAll(
    '.highlight-bash .highlight pre, .highlight-python .highlight pre'
  );
  codeBlocks.forEach((codeBlock) => {
    initializeCodeBlock(codeBlock);
  });
});
