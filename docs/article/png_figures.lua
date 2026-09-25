-- Word n'affiche pas les figures PDF : on leur substitue les PNG de figures/.
function Image(img)
  img.src = "../../../figures/" .. img.src:gsub("%.pdf$", ".png")
  return img
end
