import argparse
from Bio import SeqIO
from PIL import Image, ImageDraw, ImageFont


Stops = set(["TAA", "TAG", "TGA"])
Colors = ["Red", "Green", "Blue"]

Starts = set(["ATG"])
StartColor = ["peru"]
GapColor = "Gray"
ORFColor = "steelblue"

DataHolder = []

###############################
# Default font. 
# OpenSans-Regular.ttf is used for the paper
FontFile = "Pillow/Tests/fonts/DejaVuSans.ttf"

YLabelFontSize    = 32
YLabelFontSmallSize = 21
HeaderFontSize    = 50
ProteinFontSize    = 25

PaddingLeft = 250
PaddingRight = 100

PaddingTop = 100
PaddingBottom = 100
PaddingBetween = 80

PanelOverhang = 50
TickSize  = 10
DotSize = 3

DpiUsed = 300

ORFWidth = 20


## Enter ORF in the format of:
##  FRAME_IN_ALN..START...END...TEXT
## Example: "0..1..100..Example"
PredictedOrfPositions = ["0..160..849..NorwaORF", "0..862..2520..Nucleoprotein"]

##############################################
# Fetching data
# Return sequence number
def _seqUploader(name):
    count = 0
    for sequence in SeqIO.parse(name, "fasta"):
        count += 1
        DataHolder.append(str(sequence.seq))
    return count

def _drawSingleStop(draw_on, x, y, size, frame):
    circle_xy   = (x, y, x + size, y + size)
    draw_on.ellipse(circle_xy, Colors[frame], outline='Black', width=0)
    return

def _drawSingleStart(draw_on, x, y, size):
    circle_xy   = (x, y, x + size, y + size)
    draw_on.ellipse(circle_xy, StartColor[0], outline='Black', width=0)
    return

def _drawSingleGap(draw_on, x, y, size):
    circle_xy   = (x, y, x + size, y + size)
    draw_on.ellipse(circle_xy, GapColor, outline='Black', width=0)
    return

def _drawRotatedText(canvas, angle, coord, text, font):
    img_txt = Image.new("RGBA", font.getsize(text), (255,255,255,255))
    draw_txt = ImageDraw.Draw(img_txt)
    draw_txt.text((0,0), text, font=font, fill="Black")
    w=img_txt.rotate(angle,  expand=1)
    canvas.paste(w, coord,  w)
    return

def _orfYCoord(frame):
    return PaddingTop + frame * SequenceCount * DotSize + PaddingBetween * (frame + 0.5)


## Drawing suggested ORF
## Draw according to alignment coordinates provided.
## Beware that no checkes is performed during drawing, ie it will draw even if ORF contains stops etc
## It is just for illustration purposes
## supports only forward ORFs (ie same direction as alignment)
def _drawORFs(drawObj):
    containsorf = [False, False, False]
    for entry in PredictedOrfPositions:
        values = entry.split("..")
        y_coord = _orfYCoord(int(values[0]))
        x1_coord = PaddingLeft + int(values[1]) * DotSize / 3
        x2_coord = PaddingLeft + int(values[2]) * DotSize / 3
        x2_line = x2_coord - ORFWidth * 3
        x2_rect = x2_coord - ORFWidth * 2  
        drawObj.rectangle([(x1_coord, y_coord - ORFWidth), (x2_rect, y_coord + ORFWidth)], ORFColor, outline="Black")
        drawObj.regular_polygon(bounding_circle=(x2_rect, y_coord, ORFWidth * 2), n_sides=3, rotation=270, fill=ORFColor, outline="Black")
        drawObj.line([(x2_line, y_coord - ORFWidth), (x2_line, y_coord + ORFWidth)], ORFColor, width=1)
        # отрисовка подписей
        drawObj.text((int((x1_coord + x2_coord) / 2), y_coord),  values[3], fill="White", font=ProteinFont, anchor="mm")
        containsorf[int(values[0])] = True
    for i in range(0, 3):
        if (not containsorf[i]):
            y_coord = _orfYCoord(i) - 25
            Vizualization.text((TripletsAmount * DotSize / 2, y_coord), "No ORFs predicted in Frame +" + str(i + 1), fill="Black", font=YLabelFont, anchor="la")

    return


##############################################
# Parameters
parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", type=str, help="Input file. Expected to be a proper alignment", required=True)
parser.add_argument("-o", "--output", type=str, help="Name of the output file. Expected to end in .png")
parser.add_argument("-l", "--label", type=str, default="Visualization of the stop codons in the alignment", help="Image label")
parser.add_argument("-font", "--font", type=str, help="Custom font file to use")
parser.add_argument("-stop", "--stops", action="store_true", help="Visualize stops in alignment", default=False)
parser.add_argument("-start", "--starts", action="store_true", help="Visualize starts in alignment", default=False)
parser.add_argument("-gap", "--gaps", action="store_true", help="Visualize gaps in alignment", default=False)

FileName = parser.parse_args().input
if (parser.parse_args().font != None): FontFile = parser.parse_args().font

if (parser.parse_args().output != None):
    OutFileName = parser.parse_args().output
else:
    OutFileName = FileName + ".output.tiff"

ImageLabel = parser.parse_args().label


###############################################
# What to draw. If nothing is provided, draw stops
DrawStops = parser.parse_args().stops
DrawStarts = parser.parse_args().starts
DrawGaps = parser.parse_args().gaps
if (DrawStops):
    YRotLabel = "In-frame stops in each sequence"
elif(DrawGaps):
    YRotLabel = "In-frame gaps in each sequence"
elif(DrawStarts):
    YRotLabel = "In-frame starts in each sequence"
if (not DrawGaps and not DrawStarts and not DrawStops):
    DrawStops = True
    YRotLabel = "In-frame stops in each sequence"

SequenceCount = _seqUploader(FileName)
print(SequenceCount)
SequenceCount += 2 # + 2 for drawing margin
TripletsAmount = int(len(DataHolder[0]) / 3)


#############################################
# Image creation

XSize = PaddingLeft + PaddingRight + TripletsAmount * DotSize
YSize = PaddingTop + PaddingBottom + SequenceCount * DotSize * 3 + PaddingBetween * 3

XStartCoordinate = PaddingLeft - DotSize
YStartCoordinate = PaddingTop

XEndCoordinate = XSize - PaddingRight + DotSize
YEndCoordinate = YSize - PaddingBottom

YLabelFont    = ImageFont.truetype(FontFile, YLabelFontSize)
YLabelFontSmall = ImageFont.truetype(FontFile, YLabelFontSmallSize)
HeaderFont    = ImageFont.truetype(FontFile, HeaderFontSize)
ProteinFont    = ImageFont.truetype(FontFile, ProteinFontSize)

Canvas = Image.new("RGBA", (XSize, YSize), (255,255,255,255))
Vizualization = ImageDraw.Draw(Canvas)

############################################
# Draw some borders etc

# Vertical lines
Vizualization.line((XStartCoordinate, YStartCoordinate, XStartCoordinate, YEndCoordinate), "Black", DotSize)
Vizualization.line((XEndCoordinate, YStartCoordinate, XEndCoordinate, YEndCoordinate), "Black", DotSize)

# Horizontal lines
for i in range(0, 4):
    x_start = XStartCoordinate - PanelOverhang
    y_coord = YStartCoordinate + SequenceCount * DotSize * i + PaddingBetween * i
    Vizualization.line((x_start, y_coord, XEndCoordinate, y_coord), "Black", DotSize)
    if (i < 3):
        y_coord += PaddingBetween
        x_start += PanelOverhang
        Vizualization.line((x_start, y_coord, XEndCoordinate, y_coord), "Black", int(DotSize / 2))

############################################
# Draw alignment coordinates 
for i in range(0, len(DataHolder[0]), 500):
    x_coord = x_coord = PaddingLeft + i
    y1_coord = YEndCoordinate
    y2_coord = YEndCoordinate + TickSize
    Vizualization.line((x_coord, y1_coord, x_coord, y2_coord), "Black", DotSize * 2)
    Vizualization.text((x_coord, y2_coord + 5), str(i), fill="Black", font=YLabelFont, anchor="la")

Vizualization.text((TripletsAmount * DotSize / 2, y2_coord + 35), "Alignment coordinate (nt)", fill="Black", font=YLabelFont, anchor="la")

###############################################
## Draw Y labels
for i in range(0, 3):
    x_coord = 40
    rot_x_coord = XStartCoordinate - 40
    y_coord = int(YStartCoordinate + PaddingBetween * i + SequenceCount * DotSize * (i * 2 + 1) / 2)
    y_rot_coord1 = y_coord - 80
    y_rot_coord2 = YStartCoordinate + SequenceCount * DotSize * i + PaddingBetween * i + 10
    Vizualization.text((x_coord, y_coord), "Frame +" + str(i + 1), fill="Black", font=YLabelFont, anchor="la")
    _drawRotatedText(Canvas, 90, (rot_x_coord, y_rot_coord1), YRotLabel, YLabelFontSmall)
    _drawRotatedText(Canvas, 90, (rot_x_coord, y_rot_coord2), "ORFs", YLabelFontSmall)

Vizualization.text((XSize / 2 - 500, 15), ImageLabel, fill="Black", font=HeaderFont, anchor="la")



############################################
# Draw ORFs
_drawORFs(Vizualization)

#############################################
# Stop codon determination and drawing
#
# ToDo: Handling of the ambiguous NT? Ask GG
number = -1
for seq in DataHolder:
    number += 1 
    for j in range (0, 3):
        y_coord = PaddingTop + j * SequenceCount * DotSize + number * DotSize + PaddingBetween * (j + 1)
        num = -1    
        for i in range(0, len(seq) - 3, 3):
            num += 1
            triplet = seq[i + j : i + j + 3]
            x_coord = PaddingLeft + num * DotSize
            if (DrawGaps and "-" in triplet):
                _drawSingleGap(Vizualization, x_coord, y_coord, DotSize)
            if (DrawStops and triplet in Stops):
                _drawSingleStop(Vizualization, x_coord, y_coord, DotSize, j)
            if (DrawStarts and triplet in Starts):
                _drawSingleStart(Vizualization, x_coord, y_coord, DotSize)


# в конце не забыть удалить drawObject
del Vizualization
Canvas.save(OutFileName, "TIFF", dpi=(DpiUsed,DpiUsed))
