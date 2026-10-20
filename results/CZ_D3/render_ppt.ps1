# Document rendering only. Requires Microsoft PowerPoint installed locally.
$ErrorActionPreference='Stop'
$out=$PSScriptRoot
$pp=New-Object -ComObject PowerPoint.Application
try {
 $deck=$pp.Presentations.Open((Join-Path $out 'CZ_D3_INTERNAL_BRIEFING_DRAFT.pptx'),$true,$false,$false)
 $deck.Export((Join-Path $out 'slide_previews'),'PNG',1600,900)
 $deck.SaveAs((Join-Path $out 'CZ_D3_INTERNAL_BRIEFING_DRAFT.pdf'),32)
 $findings=@()
 foreach($sl in $deck.Slides) { foreach($shape in $sl.Shapes) {
  if($shape.HasTextFrame -eq -1 -and $shape.TextFrame.HasText -eq -1) {
   if($shape.TextFrame.TextRange.BoundHeight -gt ($shape.Height+2)) { $findings+=@{slide=$sl.SlideIndex;name=$shape.Name} }
  }
 } }
 @{scope='PowerPoint rendered text bounds only';findings=$findings;count=$findings.Count;not_scientific_audit=$true} | ConvertTo-Json -Depth 5 | Set-Content -Encoding UTF8 (Join-Path $out 'OFFICE_TEXT_OVERFLOW_CHECK.json')
 $deck.Close()
 [Runtime.InteropServices.Marshal]::ReleaseComObject($deck) | Out-Null
} finally { $pp.Quit(); [Runtime.InteropServices.Marshal]::ReleaseComObject($pp) | Out-Null }