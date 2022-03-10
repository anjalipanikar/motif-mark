import cairo
import re
import argparse
import matplotlib.colors as mplc


def get_args():                                                      #function to parse command line arguments
    parser = argparse.ArgumentParser(description=" ")
    parser.add_argument("-f", "--fastafile", help="input fasta file", required=True)
    parser.add_argument("-m", "--motiffile", help="input motif file", required=True)
    return parser.parse_args() 


def get_motifs(motif_file):                                          # function to retrieve motifs from file and translate to regex
    
    motifs = [] 
    
    with open(motif_file, 'r') as mh:
                                                                    
        for line in mh:
            motifs.append(line.strip())

    # list of motif matches
    motif_regex_list = []

    # nucleotide codes
    codes = {'y': '[ctu]',
             'r': '[ag]',
             's': '[cg]',
             'w': '[atu]',
             'k': '[gtu]',
             'm': '[ca]',
             'b': '[cgtu]',
             'd': '[agtu]',
             'v': '[agc]',
             'h': '[actu]'}

    # for each motif in the list
    for motif in motifs:
        
        motif = motif.lower()

        # for each nucleotide code 
        for key in codes:
            
            ambiguous_nucleotide_match = re.finditer(r'(?i)'+key, motif)

            # replace motif with the regex match 
            for i in ambiguous_nucleotide_match:
                motif = motif.replace(i.group(), codes[key])
            
        # add to motif match list
        motif_regex_list.append(motif) 

    return motifs, motif_regex_list


def get_seq(fastafile):                                         # function to get sequences from fasta file 

    with open(fastafile, 'r') as fh:

        # dict to store seqs
        seq_dict = {}
        seq = ''
        header = ''

        for line in fh:
            line = line.strip()

            # key = header
            if line[0] == '>':
                header = line
                seq = ''

            # value = seq
            elif line[0] != '>':
                seq += line
                seq_dict[header] = seq

    return seq_dict

class AnnotateSeq:                                              # function to get sequence info 

    def __init__(self, seq):

        # sequence of interest
        self.seq = seq

        # initialize lists 
        self.indexes = []
        self.nucleotide = []

        # loop through seq and fill lists 
        for i, nuc in enumerate(self.seq):

            # check for exons 
            if nuc.isupper():
                self.indexes.append(i)
                self.nucleotide.append(nuc)

        # return index range of exon
        self.exon_locations = (self.indexes[0], self.indexes[-1])

    def motif_coords(self, motif_regex):                        # get start position for motifs 

        self.motif_coords_list = []

        match = re.finditer(r'(?i)'+motif_regex, self.seq)
        for i in match:
            self.motif_coords_list.append(str(i.start()))

        return self.motif_coords_list

def get_coords():
    
    args = get_args() 

    # read in each fasta entry to a dictonary
    seq_dict = get_seq(args.fastafile)
    
    # get motifs and motif regex matches 
    motif_list, motif_matches = get_motifs(args.motiffile)
    
    temp = args.fastafile

    # outfile name
    outfile = temp.strip(".fasta") + ".png"
    
    colors = ['dodgerblue', 'orange', 'limegreen', 'fuchsia', 'red']

    # initialize lists
    exon_locations = []
    motif_coords = []
    seq_lengths = []
    
    for seq_entry in seq_dict:

        # seq obj for every fasta entry 
        seq_obj = AnnotateSeq(seq_dict[seq_entry])

        # add to exon coord list 
        exon_locations.append(seq_obj.exon_locations)

        # add to seq length list 
        seq_lengths.append(len(seq_dict[seq_entry]))

        # list to hold coords of motifs for each seq 
        seq_motifs = []

        # loop over regex motif matches
        for match in motif_matches:

            # tuple to hold pos of motifs 
            seq_motifs.append(
                tuple(seq_obj.motif_coords(match)))

        # list to hold all motif coords per entry 
        motif_coords.append(seq_motifs)

    return (motif_coords, exon_locations, seq_lengths, motif_list, outfile, colors)

class pycairo():                                                        # pycairo class to create image
    def __init__(self, motifs, lengths, exons):
        self.motifs = motifs
        self.lengths = lengths
        self.exons = exons
        self.outfile = get_coords()[4]
        self.colors = get_coords()[5]

    def init_surface(self):                                             # initialize surface

        # dimensions
        n = len(self.exons)
        self.width = 100 * n
        self.height = (100 * n) + 25 * n
        self.svg_len = max(self.lengths)

        # surface
        self.surface = cairo.SVGSurface(self.outfile, 300, self.height)
        self.context = cairo.Context(self.surface)
        self.context.scale(300, 300)

        # centers
        self.centers = []
        top = 1
        self.centers.append(top/(n+1) + 0.25)

        while top < n:
            top += 1
            self.centers.append(top/(n+1) + 0.25)

    def draw_sequences(self):                                           # function to draw sequence lines 
        # create lines for each sequence 
        for l, c in zip(self.lengths, self.centers):
            x, x2 = 0, l/self.svg_len
            y, y2 = c, c
            r, g, b, a = mplc.to_rgba('dimgrey')
            self.context.set_source_rgba(r, g, b, a)
            self.context.move_to(x, y)
            self.context.set_line_width(0.01)
            self.context.line_to(x2, y2)

        self.context.stroke()

    def draw_exons(self):                                               # function to draw the exons 

        for seqlen, ecenter, exon in zip(
                self.lengths, self.centers, self.exons):

            # set exon color 
            r, g, b, a = mplc.to_rgba('dimgrey')
            self.context.set_source_rgba(r, g, b, 0.7)

            # draw rectangle 
            self.context.rectangle(exon[0]/self.svg_len,
                                   (ecenter-0.02),
                                   (exon[1]-exon[0])/self.svg_len, 0.04)
            self.context.fill()

        self.context.stroke()

    def draw_motifs(self):

        colors = ['dodgerblue', 'orange', 'limegreen', 'fuchsia', 'red']

        for mot, center in zip(self.motifs, self.centers):

            # for every seq motif, length, and color 
            for m, seqlen, col in zip(mot, self.lengths, colors):

                r, g, b, a = mplc.to_rgba(col)

                # set color as context 
                self.context.set_source_rgba(r, g, b, a)

                for pos in m: #for position in motif 

                    motif = int(pos)

                    # scale to desired length 
                    x = motif/self.svg_len

                    # center the line 
                    y = center + 20/self.svg_len
                    y1 = center - 20/self.svg_len

                    # draw line
                    self.context.move_to(x, y)
                    self.context.line_to(x, y1)

                    self.context.set_line_width(0.005)

                    self.context.stroke()

        #create the legend 
        args = get_args()

        self.context.set_source_rgba(0, 0, 0, 1)
        
        self.context.select_font_face("Courier", cairo.FONT_SLANT_NORMAL,
                                      cairo.FONT_WEIGHT_BOLD)
        self.context.set_font_size(0.045)

        # x, y coordinate of text
        x = 0.024
        y = 0.1
        self.context.move_to(x, y)
        self.context.show_text("Motif Legend")
        self.context.stroke()

        # get motif text from file
        motif_txt = get_coords()[3]
        self.context.set_font_size(0.035)

        # draw each motif with its color in the rectangle
        for i, mtext in enumerate(motif_txt):
            r, b, g, a = mplc.to_rgba(colors[i])
            self.context.set_source_rgba(r, b, g, a)

            self.context.move_to(x, 0.15 + 0.05*i)
            self.context.show_text(mtext)
            self.context.stroke()
            
        seq_dict = get_seq(args.fastafile)
        
        x = 0.020
        y = 0.4 
        
        self.context.set_source_rgb(0, 0, 0)
        self.context.set_font_size(0.025)
        
        for key in seq_dict.keys():
            self.context.move_to(x, y)
            key = key.strip('>')
            self.context.show_text(key)
            y += 0.2

    def finish_surface(self):

        # start svg
        svg = pycairo(self.motifs, self.lengths, self.exons)

        # setup svg surface
        svg.init_surface()

        # draw the scaled regions for each entry
        svg.draw_sequences()

        # draw exons
        svg.draw_exons()

        svg.draw_motifs()

        # finish the surface!
        svg.surface.finish()


def main():

    #collect all necessary seq information
    motifs, exons, lengths = get_coords()[0:3]

    #draw the image 
    pycairo(motifs, lengths, exons).finish_surface()


main()
    
    

