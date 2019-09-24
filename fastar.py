###
#
# === FASTAr (FASTA Processor) ===
#
#   = Written by Jake Bradford
#   = Licensed under MIT license. See LICENSE file.
#
#   = See README.md for detailed instructions and verification method.
#   = Or use: python fastar.py --help
#
# === END ===
#
###

import argparse, os, time, string


## Globals declaration
MODES = ['collapse', 'extract', 'analyse', 'refadjust', 'refgeneextract', 'singularise']

def positionAdjust(position, offset, lineLength):
    return max(0, int(position) - offset)# - int(position) / lineLength)

def adjustStringArray(inputArray, offset, lineLength):
    return '%s,' % ','.join(map(str, [positionAdjust(x, offset, lineLength) for x in map(int, inputArray[:-1].split(','))]))
    

def main(args):
    global MODES

    ### Args validation
    if args.m not in MODES:
        print 'Invalid mode: %s\n\tChoose from: %s' % (args.m, ','.join(MODES))
        exit
    if os.path.exists(args.f) == False:
        print 'Could not find file: %s' % args.f
        exit

    ### Collapse FASTA into single line
    if args.m == 'collapse':
        print 'Collapsing into single line (but leaving FASTA comments in tact)'
        
        if not args.upper:
            print '\n\n=== COLLAPSING WITHOUT CAPITALISING ===\n\n'
            time.sleep(1)
        
        with open(args.f, 'r') as fRead:
            with open('%s.collapse' % args.f, 'w+') as fWrite:
                fRead = fRead.read()
                fRead = fRead.replace("\r", "\n").replace("\n\n", "\n")
                lines = fRead.split("\n")
                seqLine = ""
                print 'Total lines: %s' % len(lines)
                i = 0
                for line in lines:
                    line = line.strip()
                    if len(line) > 1 and line[0] == ">":  # header
                        fWrite.write('%s\n' % line)
                        if args.upper == True:
                            seqLine = seqLine.upper()
                        fWrite.write('%s' % seqLine)
                        seqLine = ""
                    else:
                        seqLine = seqLine + line
                    if i % 10000 == 0:
                        print 'Progress: %s completed (of %s)' % (i, len(lines))
                        if args.upper == True:
                            seqLine = seqLine.upper()
                        fWrite.write('%s' % seqLine)
                        fWrite.flush()
                        seqLine = ""
                    i = i + 1
                print 'Finished %s of %s' % (i, len(lines))
                if args.upper == True:
                    seqLine = seqLine.upper()        
                fWrite.write(seqLine)
        print 'Done.'

    ### Extract sub-sequence
    if args.m == 'extract':
        start = int(args.s)
        end = int(args.e)
        with open(args.f, 'r') as fRead:
            with open('%s.extract' % args.f, 'w+') as fWrite:
                fRead = fRead.read()
                fRead = fRead.replace("\r", "\n").replace("\n\n", "\n")
                lines = fRead.split("\n")
                
                print 'Removing FASTA comment lines (lines starting with >)...'
                i = 0
                headerLine = lines[0]
                while i < len(lines):
                    line = lines[i].strip()
                    if len(lines[i]) > 1 and lines[i][0] == ">":
                        print 'Deleting FASTA comment: %s' % lines[i]
                        del lines[i]
                    elif len(lines[i]) == 0:
                        print 'Deleting empty line #%s' % i
                        del lines[i]
                    i = i + 1

                print 'Collapsing file...'
                sequence = ''.join(lines).replace('\n', '').replace('\r', '')
                
                print 'Input sequence length after processing: %s' % len(sequence)
                print 'Output sequence length: %s' % (len(sequence[(start):end]))
                print 'Output starts with: %s\nOutput ends with: %s' % ((sequence[(start - 1):end])[:10], (sequence[(start - 1):end])[-10:])
                print 'Saving output...'
                fWrite.write('%s-extract[%s-%s]\n' % (headerLine, args.s, args.e))
                if args.upper:
                    sequence = sequence.upper()
                fWrite.write(sequence[(start):end])
        print 'Done.'

    ### Adjust a FASTA annotation file
    if args.m == 'refadjust':
        print args
        offset = int(args.o)
        lineLength = int(args.l)
        if args.c != None:
            restrictToChrm = args.c.split(',')
        else:
            restrictToChrm = None
        with open(args.f, 'r') as fRead:
            with open('%s.refadjust' % args.f, 'w+') as fWrite:
                fRead = fRead.read()
                fRead = fRead.replace("\r", "\n").replace("\n\n", "\n")
                lines = fRead.split("\n")
                seqLine = ""
                print 'Total lines: %s' % len(lines)
                i = 0
                # write and remove the header line
                for line in lines:       
                    data = line.split('\t')
                    if len(data) < 10:
                        continue
                    
                    # keep the original length of this gene
                    originalLength = int(data[5]) - int(data[4])
                        
                    # we want to modify txStart (4), txEnd (5), cdsStart (6), cdsEnd (7) (integers)
                    # and exonStarts(9), exonEnds (10) (arrays).
                    # if val < 0: return 0
                    data[4] = positionAdjust(data[4], offset, lineLength)
                    data[5] = positionAdjust(data[5], offset, lineLength)
                    data[6] = positionAdjust(data[6], offset, lineLength)
                    data[7] = positionAdjust(data[7], offset, lineLength)

                    # calculate new length
                    newLength = int(data[5]) - int(data[4])
                    
                    # split CSV (remove last comma), convert to INT, foreach calculate offset, convert to str, construct CSV with extra comma    
                    data[9] = adjustStringArray(data[9], offset, lineLength) #'%s,' % ','.join(map(str, [max(0, x - offset) for x in map(int, data[9][:-1].split(','))]))
                    data[10] = adjustStringArray(data[10], offset, lineLength) #'%s,' % ','.join(map(str, [max(0, x - offset) for x in map(int, data[10][:-1].split(','))]))

                    if i % 100 == 0:
                        print 'Progress: %s completed (of %s)' % (i, len(lines))
                    
                    # Determine whether the annotation line should be printed
                    keepDifferentLength = (newLength != originalLength and args.removeDifferentLength == False)
                    isSameLength = (newLength == originalLength)
                    isWithinDesiredSetOfChr = (restrictToChrm == None or data[2] in restrictToChrm)
                    isWithinBoundsOfAdjustment = (args.s == None or args.e == None) or (args.removeOutsideBounds and int(args.s) <= int(data[4]) and int(data[5]) <= int(args.e))
                    
                    if (keepDifferentLength or isSameLength) and isWithinDesiredSetOfChr and isWithinBoundsOfAdjustment:
                        fWrite.write('%s\n' % ('\t'.join(map(str, data))))
                    i = i + 1
                
                print 'Finished %s of %s' % (i, len(lines))
                
        print 'Done.'                    

    ### Adjust a FASTA annotation file
    if args.m == 'analyse':
        with open(args.f, 'r') as fRead:
            print '\nAnalysing file: %s' % args.f
            fRead = fRead.read()
            print '\nTotal file length: %s' % len(fRead)
            fRead = fRead.replace("\r", "\n").replace("\n\n", "\n")
            lines = fRead.split("\n")
            print 'Total lines in file: %s' % len(lines)
            seqLine = ""
            i = 0
            fastaBlockLines = 0
            fastaBlockLength = 0
            lastFastaComment = ""
            totalLength = 0
            first = True
            for line in lines: 
                if len(line) == 0:
                    print '\nEmpty line'
                    continue
                if len(line) > 0 and line[0] == ">":
                    if not first:
                        print '\nFasta block beginning: %s' % lastFastaComment
                        print 'Block length: %s chars' % fastaBlockLength
                        print 'Block length: %s lines' % fastaBlockLines
                    lastFastaComment = line
                    fastaBlockLines = 0
                    fastaBlockLength = 0
                else:
                    first = False
                    fastaBlockLines = fastaBlockLines + 1
                    fastaBlockLength = fastaBlockLength + len(line)
                    totalLength = totalLength + len(line)
                    
            print '\nFasta block beginning: %s' % lastFastaComment
            print 'Block length: %s chars' % fastaBlockLength
            print 'Block length: %s lines' % fastaBlockLines
            
            print '\nTotal file length (excl. FA comments): %s' % totalLength
            
            print 'Done.'
    
    ### Extract gene names between start and end position using the given reference
    if args.m == 'refgeneextract':
        start = int(args.s)
        end = int(args.e)

        with open(args.f, 'r') as fRead:
            with open('%s.genes' % args.f, 'w+') as fWrite:
                for line in fRead.split('\n'):
                    lineSplit = line.split('\t')
                    if int(lineSplit[4]) >= start and int(lineSplit[5]) <= end:
                        fWrite.write('%s\n' % lineSplit[8])
    
    ### Convert a multi-fasta file to multiple single fasta files
    if args.m == 'singularise':
        with open(args.f, 'r') as fRead:
            fWrite = None
            fWriteDir = os.path.dirname(args.f)
            
            # for each line in the file
            for line in fRead:
                line = line.strip()
                
                # just found a new fasta segment. open a new file
                if line[0] == '>':
                
                    # prepare the file name
                    charsToSwitch = ' +-.:;'
                    processedName = line[1:].translate(string.maketrans(charsToSwitch, '_' * len(charsToSwitch)))
                    fileName = os.path.join(fWriteDir, '%s.fa' % processedName)
                    
                    # print the new file name as a progress indicator
                    print fileName
                    
                    # close the current file if necessary
                    if fWrite is not None:
                        fWrite.close()
                        
                    # open a new one
                    fWrite = open(fileName, 'w+')

                # write the line to the file
                if line[0] != '>':
                    fWrite.write('%s\n' % line)
            
            # close the last file if necessary
            if fWrite is not None:
                fWrite.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = 'FASTAr (FASTA Processor)')
    parser.add_argument('-m', help='Mode (%s)' % ' | '.join(MODES), default=None, required=True)
    parser.add_argument('-f', help='File to process', default=None, required=True)
    parser.add_argument('-s', help='[extract/refgeneextract] start position (inclusive)', default=None)
    parser.add_argument('-e', help='[extract/refgeneextract] end position (inclusive)', default=None)
    parser.add_argument('-o', help='[refadjust] offset (x - offset)', default=None)
    parser.add_argument('-l', help='[refadjust] line length of original file (UCSC = 50) (does nothing)', default=50)
    parser.add_argument('-c', help='[refadjust] include only these chromosomes as CSV format (eg: chr1,chr2,chr3) (default: all)', default=None)
    parser.add_argument('--removeDifferentLength', help='[refadjust] Remove annotation if transcription LENGTH differs after adjustment. It would change if the transcription is trimmed because it was outside the bounds of -s and -e (after -o).', const=True, default=False, nargs='?')
    parser.add_argument('--removeOutsideBounds', help='[refadjust] Remove annotation if it would run outside the bounds of -s to -e.', const=True, default=False, nargs='?')
    parser.add_argument('--upper', help='[collapse] Convert output sequence to uppercase', const=True, default=False, nargs='?')
    args = parser.parse_args()
    
    if os.path.exists(args.f) == False:
        print 'File provided does not exist: %s' % args.f
        exit()
    
    if args.m in ['extract', 'refgeneextract']:
        if (args.s == None or args.e == None):
            print 'Start (%s) and end (%s) numbers cannot be None.' %(args.s, args.e)
            exit()
        
        if (int(args.s) >= int(args.e)):
            print 'End number must be greater than start number.'
            exit()
    
    main(args)
