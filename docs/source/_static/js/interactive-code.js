document.addEventListener('DOMContentLoaded', () => {
  const codeBlock = document.querySelector('.highlight-bash .highlight pre');
  if (codeBlock) {
    codeBlock.contentEditable = true;
    codeBlock.spellcheck = false;
    codeBlock.style.outline = 'none';
    codeBlock.style.color = 'black';
    codeBlock.style.cursor = 'text';
    codeBlock.style.userSelect = 'text';
    codeBlock.style.whiteSpace = 'pre';
    codeBlock.style.paddingTop = '24px';

    // Sanitize input to plain text only
    codeBlock.addEventListener('paste', (e) => {
      e.preventDefault();
      const text = e.clipboardData.getData('text/plain');
      document.execCommand('insertText', false, text);
    });

    // Prevent HTML insertion
    codeBlock.addEventListener('input', (e) => {
      if (e.inputType === 'insertFromPaste') {
        // Already handled by paste event
        return;
      }
      // For other input, ensure only text content
      const textContent = codeBlock.textContent;
      if (codeBlock.innerHTML !== textContent) {
        codeBlock.textContent = textContent;
      }
    });

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
      color: #666;
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
      color: #666;
    `;
    parent.appendChild(resetButton);

    // Add bash tag
    const bashTag = document.createElement('div');
    bashTag.textContent = 'bash';
    bashTag.style.cssText = `
      position: absolute;
      top: 0;
      left: 0;
      background: transparent;
      color: #5f6368;
      padding: 4px 8px;
      font-size: 11px;
      font-weight: 500;
      z-index: 5;
    `;
    parent.appendChild(bashTag);

    codeBlock.addEventListener('click', () => {
      codeBlock.style.border = '1px solid black';
    });

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

    // Store original text
    const originalText = codeBlock.textContent;

    // Copy functionality
    copyButton.addEventListener('click', (e) => {
      e.stopPropagation();
      navigator.clipboard.writeText(codeBlock.textContent);

      // Show message and keep it visible
      message.style.opacity = '1';
    });

    // Reset functionality
    resetButton.addEventListener('click', (e) => {
      e.stopPropagation();
      codeBlock.textContent = originalText;
      message.style.opacity = '0';
    });
  }
});
