import numpy as np
import copy

def search_aa(seqlen, peptide_sequence, searched_seqs, locked_positions=None):
    '''
    Search for a new amino acid sequence
    '''
    restypes = np.array(['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I',
                 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V'])
    #Go through a shuffled version of the positions and aas
    for pi in np.random.choice(np.arange(seqlen), seqlen, replace=False):
        if locked_positions and pi in locked_positions:
            continue  # Skip locked positions
        for aa in np.random.choice(restypes, len(restypes), replace=False):
            new_seq = list(peptide_sequence)
            new_seq[pi] = aa
            new_seq = ''.join(new_seq)
            if new_seq in searched_seqs:
                #Recursion - next will use new seq as peptide_sequence
                return search_aa(seqlen, new_seq, searched_seqs, locked_positions)
            else:
                # Ensure locked positions are still correct
                if locked_positions:
                    for pos, locked_aa in locked_positions.items():
                        if new_seq[pos] != locked_aa:
                            break
                    else:
                        return new_seq
                else:
                    return new_seq
    return None  # Return None if no valid sequence is found

def mutate_sequence(peptide_sequence, sequence_scores, locked_positions=None):
    '''Mutate the amino acid sequence randomly, respecting locked positions
    '''
    restypes = np.array(['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I',
                 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V'])
    seqlen = len(peptide_sequence)
    searched_seqs = sequence_scores['sequence']
    
    #Mutate seq
    seeds = [peptide_sequence]
    #Go through a shuffled version of the positions and aas
    for seed in seeds:
        for pi in np.random.choice(np.arange(seqlen), seqlen, replace=False):
            if locked_positions and pi in locked_positions:
                continue  # Skip locked positions
            for aa in np.random.choice(restypes, len(restypes), replace=False):
                new_seq = list(seed)
                new_seq[pi] = aa
                new_seq = ''.join(new_seq)
                if new_seq in searched_seqs:
                    continue
                else:
                    # Ensure locked positions are still correct
                    if locked_positions:
                        for pos, locked_aa in locked_positions.items():
                            if new_seq[pos] != locked_aa:
                                break
                        else:
                            return new_seq
                    else:
                        return new_seq
        seeds.append(new_seq)
    
    # If no valid mutation is found, try searching for a completely new sequence
    return search_aa(seqlen, peptide_sequence, searched_seqs, locked_positions)
