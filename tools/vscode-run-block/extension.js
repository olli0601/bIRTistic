// Minimal VS Code extension: run the block bounded by '# ----' (or '# %%')
// markers around the cursor, in the Python terminal, via the Python extension's
// execSelectionInTerminal. Bound to Cmd+Shift+Enter (see package.json).
const vscode = require('vscode');

const MARKER = /^\s*#\s*(----+|%%)/;   // a block boundary line

function activate(context) {
  const cmd = vscode.commands.registerCommand('runBlockToMarker.run', async () => {
    const ed = vscode.window.activeTextEditor;
    if (!ed) return;
    const doc = ed.document;
    const cur = ed.selection.active.line;
    const onMarker = MARKER.test(doc.lineAt(cur).text);

    // start = line after the marker at/above the cursor (else top of file)
    let start = 0;
    for (let i = onMarker ? cur : cur - 1; i >= 0; i--) {
      if (MARKER.test(doc.lineAt(i).text)) { start = i + 1; break; }
    }
    // end = line before the next marker below (else end of file)
    let end = doc.lineCount - 1;
    for (let i = cur + 1; i < doc.lineCount; i++) {
      if (MARKER.test(doc.lineAt(i).text)) { end = i - 1; break; }
    }
    // trim trailing blank lines
    while (end > start && doc.lineAt(end).text.trim() === '') end--;
    if (end < start) return;

    ed.selection = new vscode.Selection(
      new vscode.Position(start, 0),
      new vscode.Position(end, doc.lineAt(end).text.length),
    );
    await vscode.commands.executeCommand('python.execSelectionInTerminal');
  });
  context.subscriptions.push(cmd);
}

function deactivate() {}

module.exports = { activate, deactivate };
