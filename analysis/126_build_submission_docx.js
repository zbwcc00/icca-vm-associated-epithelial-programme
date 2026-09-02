const fs = require('fs');
const path = require('path');
const {
  Document, Packer, Paragraph, TextRun, HeadingLevel, PageBreak,
  ImageRun, AlignmentType, Footer, PageNumber, NumberFormat
} = require('docx');

const project = process.env.ICCA_VM_PROJECT_ROOT || path.resolve(__dirname, '..');
const markdownPath = path.join(project, 'manuscript_draft', '09_full_manuscript_draft_v1.md');
const outputPath = path.join(project, 'manuscript_draft', '09_full_manuscript_submission_v1.docx');
const figureDir = path.join(project, 'figures', 'manuscript_v1');

function cleanMarkdown(text) {
  return text
    .replace(/\*\*([^*]+)\*\*/g, '$1')
    .replace(/\*([^*]+)\*/g, '$1')
    .replace(/`([^`]+)`/g, '$1')
    .replace(/\u2014/g, '—');
}

function paragraphFromText(text, options = {}) {
  return new Paragraph({
    alignment: options.alignment,
    spacing: { after: options.after ?? 120, line: 276 },
    indent: options.indent,
    children: [new TextRun({ text: cleanMarkdown(text), font: 'Arial', size: options.size ?? 22 })]
  });
}

function pngDimensions(imageData) {
  const signature = '89504e470d0a1a0a';
  if (imageData.subarray(0, 8).toString('hex') !== signature) {
    throw new Error('Expected a PNG figure asset');
  }
  return { width: imageData.readUInt32BE(16), height: imageData.readUInt32BE(20) };
}

function fittedImageSize(imageData, maxWidth = 900, maxHeight = 650) {
  const { width, height } = pngDimensions(imageData);
  const scale = Math.min(maxWidth / width, maxHeight / height);
  return { width: Math.round(width * scale), height: Math.round(height * scale) };
}

const markdown = fs.readFileSync(markdownPath, 'utf8').replace(/\r\n/g, '\n');
const lines = markdown.split('\n');
const children = [];
let inReferences = false;

for (const line of lines) {
  if (!line.trim()) {
    children.push(new Paragraph({ spacing: { after: 80 } }));
    continue;
  }
  if (line.startsWith('# ')) {
    children.push(new Paragraph({ alignment: AlignmentType.CENTER, spacing: { after: 240 }, children: [new TextRun({ text: cleanMarkdown(line.slice(2)), bold: true, font: 'Arial', size: 30 })] }));
    continue;
  }
  if (line.startsWith('## ')) {
    const heading = line.slice(3);
    if (heading === 'References') inReferences = true;
    children.push(new Paragraph({ heading: HeadingLevel.HEADING_1, spacing: { before: 240, after: 120 }, children: [new TextRun({ text: cleanMarkdown(heading), bold: true, font: 'Arial', size: 26 })] }));
    continue;
  }
  if (line.startsWith('### ')) {
    children.push(new Paragraph({ heading: HeadingLevel.HEADING_2, spacing: { before: 180, after: 100 }, children: [new TextRun({ text: cleanMarkdown(line.slice(4)), bold: true, font: 'Arial', size: 24 })] }));
    continue;
  }
  if (line.startsWith('|')) {
    children.push(paragraphFromText(line.replace(/\|/g, '  '), { size: 18, after: 60 }));
    continue;
  }
  if (/^\d+\.\s/.test(line)) {
    children.push(paragraphFromText(line, { indent: { left: 360, hanging: 240 }, size: 18, after: 70 }));
    continue;
  }
  children.push(paragraphFromText(line, { size: inReferences ? 18 : 22, after: inReferences ? 80 : 130 }));
}

const figureFiles = [
  ['Figure 1. Study design and evidence boundary', 'Fig1_study_design_and_claim_boundaries.png'],
  ['Figure 2. Cross-cohort definition and patient coverage', 'Fig2_patient_level_VMlike_coverage.png'],
  ['Figure 3. Constrained receiver-side associations', 'Fig3_constrained_receiver_associations.png'],
  ['Figure 4. Directional replication of bulk associations', 'Fig4_three_bulk_directional_replication.png'],
  ['Figure 5. Matched RNA–protein concordance', 'Fig5_matched_proteomic_concordance.png'],
  ['Figure 6. Association-only genomic context', 'Fig6_association_only_genomic_context.png']
];

children.push(new Paragraph({ children: [new PageBreak()] }));
children.push(new Paragraph({ heading: HeadingLevel.HEADING_1, children: [new TextRun({ text: 'Figures', bold: true, font: 'Arial', size: 28 })] }));
for (const [title, filename] of figureFiles) {
  children.push(new Paragraph({ heading: HeadingLevel.HEADING_2, spacing: { before: 240, after: 120 }, children: [new TextRun({ text: title, bold: true, font: 'Arial', size: 24 })] }));
  const imagePath = path.join(figureDir, filename);
  if (fs.existsSync(imagePath)) {
    const imageData = fs.readFileSync(imagePath);
    children.push(new Paragraph({ alignment: AlignmentType.CENTER, spacing: { after: 160 }, children: [new ImageRun({ type: 'png', data: imageData, transformation: fittedImageSize(imageData) })] }));
  } else {
    children.push(paragraphFromText(`[Figure file not found: ${filename}]`, { size: 18 }));
  }
}

const document = new Document({
  creator: 'Authors',
  title: 'VM-associated epithelial transcriptional programme in iCCA',
  description: 'Submission-mode manuscript draft',
  styles: {
    default: { document: { run: { font: 'Arial', size: 22 }, paragraph: { spacing: { line: 276, after: 120 } } } }
  },
  sections: [{
    properties: { page: { margin: { top: 1080, right: 1080, bottom: 1080, left: 1080 } } },
    footers: { default: new Footer({ children: [new Paragraph({ alignment: AlignmentType.CENTER, children: [new TextRun({ text: 'Page ', font: 'Arial', size: 18 }), new TextRun({ children: [PageNumber.CURRENT], font: 'Arial', size: 18 })] })] }) },
    children
  }]
});

Packer.toBuffer(document).then(buffer => {
  fs.writeFileSync(outputPath, buffer);
  console.log(outputPath);
});
