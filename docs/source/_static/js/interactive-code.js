document.addEventListener('DOMContentLoaded', () => {
  // Find the first code block in the site module
  const codeBlocks = document.querySelectorAll(
    '.highlight-bash .highlight pre'
  );

  if (codeBlocks.length > 0) {
    const firstCodeBlock = codeBlocks[0];

    // Make the code block selectable and editable
    firstCodeBlock.contentEditable = true;
    firstCodeBlock.spellcheck = false;
    firstCodeBlock.style.outline = 'none';
    firstCodeBlock.style.border = '1px solid transparent';
    firstCodeBlock.style.borderRadius = '4px';
    firstCodeBlock.style.padding = '8px';
    firstCodeBlock.style.transition = 'border-color 0.2s ease';
    firstCodeBlock.style.cursor = 'text';

    // Add focus styles and select all text on click
    firstCodeBlock.addEventListener('click', () => {
      firstCodeBlock.style.borderColor = '#007bff';
      firstCodeBlock.style.backgroundColor = '#f8f9fa';

      // Select all text when clicked
      const range = document.createRange();
      range.selectNodeContents(firstCodeBlock);
      const selection = window.getSelection();
      selection.removeAllRanges();
      selection.addRange(range);
    });

    firstCodeBlock.addEventListener('blur', () => {
      firstCodeBlock.style.borderColor = 'transparent';
      firstCodeBlock.style.backgroundColor = 'transparent';
    });

    // Add a subtle hint
    firstCodeBlock.setAttribute(
      'title',
      'Click to select and edit this command'
    );
  }
});
