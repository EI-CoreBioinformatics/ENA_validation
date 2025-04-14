#!/ei/software/cb/python_miniconda/24.9.2-0_py3.12_swarbred/x86_64/bin/python3
from collections import defaultdict
import argparse
import os

def input_gff_file(gff3):
    """
    Parse the GFF3 file into structured dictionaries for processing.
    """
    EXON = defaultdict(list)
    UTR = defaultdict(list)
    GFF3_DICT = defaultdict(list)
    HEADER = []
    transcript_representative = {}

    with open(gff3, 'r') as file:
        for line in file:
            if line.startswith('###'):
                pass
            elif line.startswith('##'):
                HEADER.append(line)
            else:
                fields = line.strip().split('\t')
                if len(fields) >= 9:
                    feature_type = fields[2]
                    attrs = dict(attr.split('=') for attr in fields[8].split(';') if '=' in attr)

                    # Parse gene ID and transcript ID
                    if feature_type in ['gene', 'ncRNA_gene']:
                        gene_id = attrs.get('ID', '')
                        transcript_id = 'NA'

                    elif feature_type in ['mRNA', 'ncRNA']:
                        transcript_id = attrs.get('ID', '')
                        note_attr = attrs.get('Note', '')
                        rep_tag = 'rep:True' in note_attr
                        transcript_representative[transcript_id] = (
                            attrs.get('representative', 'False').lower() == 'true' or rep_tag
                        )
                    else:
                        transcript_id = attrs.get('Parent', '').split(',')[0]

                    ids = attrs.get('ID', '')
                    start = int(fields[3])
                    end = int(fields[4])
                    seqid = fields[0]
                    source = fields[1]
                    score = fields[5]
                    strand = fields[6]
                    phase = fields[7]
                    attributes = fields[8]

                    representative = transcript_representative.get(transcript_id, False)

                    info = {
                        'seqid': seqid,
                        'source': source,
                        'type': feature_type,
                        'start': start,
                        'end': end,
                        'score': score,
                        'strand': strand,
                        'phase': phase,
                        'attributes': attributes,
                        'id': ids,
                        'transcript': transcript_id,
                        'representative': representative
                    }

                    if feature_type in ['three_prime_UTR', 'five_prime_UTR']:
                        UTR[gene_id].append(info)
                        GFF3_DICT[gene_id].append(info)
                    elif feature_type == 'exon':
                        EXON[gene_id].append(info)
                        GFF3_DICT[gene_id].append(info)
                    else:
                        GFF3_DICT[gene_id].append(info)

    dictionaries = {"utr": UTR, "gff3_dict": GFF3_DICT, "exon": EXON, "header": HEADER}
    print(f"Input parsed")
    return dictionaries

def get_duplicate_features(dictionary, verbose):
    """
    Identify duplicate features only between different genes.
    """
    dup_list = []
    feature_groups = defaultdict(list)

    for gene_id, features in dictionary.items():
        for feature in features:
            if feature["type"] in ["CDS", "exon", "three_prime_UTR", "five_prime_UTR", "gene", "mRNA", "ncRNA"]:
                feature_key = (feature["seqid"], feature["start"], feature["end"], feature["strand"], feature["type"])
                feature_groups[feature_key].append(gene_id)

    if verbose:
        print("\nChecking for duplicates...")
    for feature_key, gene_ids in feature_groups.items():
        unique_gene_ids = set(gene_ids)
        if len(unique_gene_ids) > 1:
            if verbose:
                print(f"Duplicate Feature Found: {feature_key} in Genes: {', '.join(unique_gene_ids)}")
            for gene_id in unique_gene_ids:
                dup_list.append(gene_id)

    if dup_list:
        if verbose:
             print(f"\nDuplicated features detected in the following genes: {', '.join(set(dup_list))}")
    else:
        if verbose:
            print("\nNo duplicate features found.")

    return list(set(dup_list))

def resolve_duplicates_with_cds_based_adjustment(dictionary, duplicate_genes, verbose):
    """
    Resolve duplicate features using CDS-based recalculations, but only for duplicate genes.
    Retain existing CDS features linked to recalculated exons.
    """
    edited_dict = defaultdict(list)
    updated_features_log = []

    for gene_id in dictionary:
        features = dictionary[gene_id]

        if gene_id in duplicate_genes:
            cds_features = [f for f in features if f['type'] == 'CDS']
            if cds_features:
                recalculated_exons = []
                recalculated_cds = []
                for i, cds in enumerate(sorted(cds_features, key=lambda x: x['start']), start=1):
                    exon_id = f"{gene_id}.recalc_exon{len(recalculated_exons) + 1}"
                    recalculated_exon = {
                        'seqid': cds['seqid'],
                        'source': cds['source'],
                        'type': 'exon',
                        'start': cds['start'],
                        'end': cds['end'],
                        'score': cds['score'],
                        'strand': cds['strand'],
                        'phase': '.',
                        'attributes': f"ID={exon_id};Parent={cds['transcript']}"
                    }
                    recalculated_exons.append(recalculated_exon)

                    recalculated_cds.append({
                        'seqid': cds['seqid'],
                        'source': cds['source'],
                        'type': 'CDS',
                        'start': cds['start'],
                        'end': cds['end'],
                        'score': cds['score'],
                        'strand': cds['strand'],
                        'phase': cds['phase'],
                        'attributes': f"ID={cds['transcript']}.CDS{i};Parent={cds['transcript']}"
                    })

                updated_features_log.extend(recalculated_exons + recalculated_cds)
                edited_dict[gene_id].extend(recalculated_exons)
                edited_dict[gene_id].extend(recalculated_cds)

                # Adjust gene and ncRNA/mRNA coordinates CHANGE
                child_features = recalculated_exons + recalculated_cds
                for feature in features:
                    if feature['type'] in ['gene', 'ncRNA']: # check david how to handle ncRNA
                        #for each gene, mrna and ncrna make a new start and end feature depending on those from the recalculated
                        feature['start'] = min(f['start'] for f in child_features)
                        feature['end'] = max(f['end'] for f in child_features)
                        edited_dict[gene_id].append(feature)
                    elif feature['type'] == 'mRNA':
                         # Collect child features (exons and CDS) that belong to this mRNA
                        child_exons_for_mrna = []
                        for f in child_features:
                            parent_id = None
                            for part in f['attributes'].split(';'):
                                if part.startswith("Parent="):
                                    parent_ids = part.split('=')[1]
                                    # Check if this mRNA's ID is among the parent IDs
                                    if feature['id'] in parent_ids:
                                        parent_id = feature['id']
                                        break
                            if parent_id == feature['id']:
                                child_exons_for_mrna.append(f)

                        # If child features are found, update mRNA coordinates
                        if child_exons_for_mrna:
                            feature['start'] = min(child['start'] for child in child_exons_for_mrna)
                            feature['end'] = max(child['end'] for child in child_exons_for_mrna)
    
                        edited_dict[gene_id].append(feature)
            else:
                edited_dict[gene_id].extend(features)
        else:
            edited_dict[gene_id].extend(features)


    # Enhanced logging
    if verbose:
        if updated_features_log:
            print(f"Updated features based on exon recalculations:")
            for feature in updated_features_log:
                print(f"  Updated Feature: {feature['attributes']} {feature['start']}-{feature['end']} ({feature['strand']})")
        else:
            print(f"No updated features based on exon recalculations.")

    return edited_dict

def extract_transcript_version(transcript):
    """ Extract the transcript version number (part after the dot) for sorting. """
    try:
        return int(transcript.split('.')[-1])  # Get the number after the last dot
    except ValueError:
        return 0  # Default to 0 if extraction fails

def check_abutting_exons(genes, verbose=True):
    """
    Identify abutting exons and enforce transcript removal.
    The goal is to remove the fewest possible transcripts while ensuring no abutting exons remain.
    """

    abutting_pairs = []
    transcript_abut_count = defaultdict(int)
    transcript_clusters = defaultdict(set)
    transcript_representative = {}

    # Step 1: Identify abutting exon pairs and group them into clusters
    if verbose:
        print("\nIdentifying abutting exons\n")
    for gene_id, exons in genes.items():
        for i, exon1 in enumerate(exons):
            for j in range(i + 1, len(exons)):
                exon2 = exons[j]
                if exon1["strand"] == exon2["strand"]:
                    if exons[i]['end'] == exons[j]['start'] - 1 or exons[i]['start'] == exons[j]['end'] + 1:
                        abutting_pairs.append((exon1, exon2))
                        transcript_abut_count[exon1["transcript"]] += 1
                        transcript_abut_count[exon2["transcript"]] += 1
                        transcript_clusters[exon1["transcript"]].add(exon2["transcript"])
                        transcript_clusters[exon2["transcript"]].add(exon1["transcript"])

                        if verbose:
                            print(f"Abutting exons found in gene {gene_id}:")
                            print(f"  {exon1['id']} ({exon1['transcript']}) : {exon1['start']} - {exon1['end']}")
                            print(f"  {exon2['id']} ({exon2['transcript']}) : {exon2['start']} - {exon2['end']}\n")

    # Step 2: Identify representative transcripts
    for gene_id, exons in genes.items():
        for exon in exons:
            transcript_representative[exon["transcript"]] = exon.get("representative", False)

#    if verbose:
#        print("\n### Debug: Checking which transcripts are representative ###\n")
#        for transcript, is_rep in transcript_representative.items():
#            print(f"Transcript: {transcript} | Representative: {is_rep}")

    # Step 3: Group all transcripts with abutting exons into clusters
    resolved_clusters = []
    visited = set()

#    if verbose:
#        print("\n### Debug: Resolving abutting transcript clusters ###\n")
    for transcript, related_transcripts in transcript_clusters.items():
        if transcript not in visited:
            cluster = set()
            queue = [transcript]
            while queue:
                curr = queue.pop()
                if curr not in visited:
                    visited.add(curr)
                    cluster.add(curr)
                    queue.extend(transcript_clusters[curr])
            resolved_clusters.append(cluster)
#            if verbose:
#                print(f"Cluster formed: {cluster}")

    # Step 4: Forced removal of transcripts
    models_to_remove = set()

    for cluster in resolved_clusters:
        rep_transcripts = {t for t in cluster if transcript_representative.get(t, False)}
        non_rep_transcripts = [t for t in cluster if t not in rep_transcripts]

#        if verbose:
#            print(f"\n### Debug: Processing cluster: {cluster} ###")
#            print(f"  Representative transcripts: {rep_transcripts}")

        if not rep_transcripts:
            # **Clusters without a representative transcript: original behavior**
            sorted_non_reps = sorted(
                non_rep_transcripts,
                key=lambda t: (-transcript_abut_count[t], -extract_transcript_version(t))
            )

            while sorted_non_reps:
                to_remove = sorted_non_reps.pop(0)
                models_to_remove.add(to_remove)

                # **Recalculate remaining abutting conflicts**
                remaining_set = cluster - models_to_remove
                remaining_abutting_pairs = [
                    (e1, e2) for e1, e2 in abutting_pairs if e1["transcript"] in remaining_set and e2["transcript"] in remaining_set
                ]

 #               if verbose:
 #                   print(f"  Removing transcript: {to_remove}")
 #                   print(f"  Remaining abutting pairs after removing {to_remove}: {remaining_abutting_pairs}")

                if not remaining_abutting_pairs:
                    break

                # **Update conflict counts dynamically**
                for e1, e2 in abutting_pairs:
                    if e1["transcript"] == to_remove:
                        transcript_abut_count[e2["transcript"]] -= 1
                    if e2["transcript"] == to_remove:
                        transcript_abut_count[e1["transcript"]] -= 1

                # **Re-sort remaining non-representative transcripts**
                sorted_non_reps = sorted(
                    [t for t in non_rep_transcripts if t not in models_to_remove],
                    key=lambda t: (-transcript_abut_count[t], -extract_transcript_version(t))
                )

        else:
            # **Clusters WITH a representative transcript: remove the most conflicting non-representative transcript(s)**
            if non_rep_transcripts:
                sorted_non_reps = sorted(
                    non_rep_transcripts,
                    key=lambda t: (-transcript_abut_count[t], -extract_transcript_version(t))
                )

#                if verbose:
#                    print(f"\n  Non-representative transcripts before removal: {sorted_non_reps}")

                while sorted_non_reps:
                    to_remove = sorted_non_reps.pop(0)
                    models_to_remove.add(to_remove)

                    # **Recalculate remaining abutting conflicts**
                    remaining_set = cluster - models_to_remove
                    remaining_abutting_pairs = [
                        (e1, e2) for e1, e2 in abutting_pairs if e1["transcript"] in remaining_set and e2["transcript"] in remaining_set
                    ]

 #                   if verbose:
 #                       print(f"  Removing transcript: {to_remove}")
 #                       print(f"  Remaining abutting pairs after removing {to_remove}: {remaining_abutting_pairs}")

                    if not remaining_abutting_pairs:
                        break

                    # **Update conflict counts dynamically**
                    for e1, e2 in abutting_pairs:
                        if e1["transcript"] == to_remove:
                            transcript_abut_count[e2["transcript"]] -= 1
                        if e2["transcript"] == to_remove:
                            transcript_abut_count[e1["transcript"]] -= 1

                    # **Re-sort remaining non-representative transcripts**
                    sorted_non_reps = sorted(
                        [t for t in non_rep_transcripts if t not in models_to_remove],
                        key=lambda t: (-transcript_abut_count[t], -extract_transcript_version(t))
                    )

    # Step 5: Ensure at least one transcript remains
    all_transcripts = set(transcript_clusters.keys())
    remaining_transcripts = all_transcripts - models_to_remove

#    if verbose:
#        print("\n### Debug: Ensuring at least one transcript remains ###")
#        print(f"  Remaining transcripts before final check: {remaining_transcripts}")

    if not remaining_transcripts and models_to_remove:
        # Identify representative transcripts in the removal list
        rep_transcripts_to_keep = [t for t in models_to_remove if transcript_representative.get(t, False)]
        
        if rep_transcripts_to_keep:
            to_keep = rep_transcripts_to_keep[0]  # Prioritize keeping a representative transcript
            print(f"  Representative transcript retained: {to_keep}")
        else:
            # Apply sorting: first by lowest abutting conflict count, then by lowest version number
#            print(f"  Sorted transcripts for potential retention: {sorted(models_to_remove, key=lambda t: (transcript_abut_count[t], extract_transcript_version(t)))}")
            to_keep = sorted(models_to_remove, key=lambda t: (transcript_abut_count[t], extract_transcript_version(t)))[0]

        models_to_remove.remove(to_keep)

        if verbose:
            print(f"  No transcripts left! Keeping: {to_keep}")

    elif not remaining_transcripts:
        if verbose:
            print("  No abutting exons found; no transcripts need to be removed.")

    return models_to_remove  # Ensure the function returns the final list

def remove_duplicate_mrna_structures(dictionary, verbose=False):
    """
    For each gene in the dictionary, check for mRNAs with identical exon structures.
    If duplicates are found, remove the mRNA feature and its associated child features,
    retaining only one mRNA per unique exon structure. Retention is determined by:
      - Keeping the transcript(s) with representative=True if available (choosing the one with the lowest version number)
      - Otherwise, retaining the mRNA with the lowest version number.
      
    Parameters:
      dictionary (dict): Gene-keyed dictionary of features (the final edited dictionary).
      verbose (bool): If True, log details about duplicate detection and removals.
    
    Returns:
      dict: The updated dictionary with duplicate mRNAs and associated features removed.
    """
    def get_exon_structure(mrna_id, features):
        """
        For the given mRNA id, return a tuple representing its exon structure.
        Each exon is represented as a (start, end) tuple.
        The exons are collected from features that either have a 'transcript' field equal
        to the mRNA id or have a "Parent=" attribute containing the mRNA id.
        """
        exons = []
        for f in features:
            if f['type'] == 'exon':
                # Check for association via the parsed 'transcript' key.
                if f.get('transcript') == mrna_id:
                    exons.append((f['start'], f['end']))
                # Also check the attributes if the 'transcript' key is not set.
                elif "Parent=" in f['attributes']:
                    # This will match lines like: "ID=...;Parent=<mrna_id>"
                    if f"Parent={mrna_id}" in f['attributes']:
                        exons.append((f['start'], f['end']))
        # Sort exons by start coordinate (works for both strands since start < end in GFF)
        return tuple(sorted(exons, key=lambda x: x[0]))

    # Process each gene in the dictionary
    for gene_id, features in dictionary.items():
        # Extract mRNA (or ncRNA) features from the gene.
        mrna_features = [f for f in features if f['type'] in ['mRNA', 'ncRNA']]
        
        # Build a mapping from mRNA id to its exon structure.
        mrna_to_structure = {}
        for mrna in mrna_features:
            mrna_id = mrna['id']
            structure = get_exon_structure(mrna_id, features)
            mrna_to_structure[mrna_id] = structure

        # Group mRNAs by their exon structure.
        structure_groups = {}
        for mrna_id, structure in mrna_to_structure.items():
            structure_groups.setdefault(structure, []).append(mrna_id)
        
        # Determine which transcripts to remove from groups with duplicate structures.
        removal_set = set()
        for structure, mrna_ids in structure_groups.items():
            if len(mrna_ids) > 1:
                # From these duplicates, prefer any with representative=True.
                rep_ids = [
                    mrna_id for mrna_id in mrna_ids 
                    if any(f for f in mrna_features if f['id'] == mrna_id and f.get('representative', False))
                ]
                if rep_ids:
                    # Retain the representative with the lowest version number.
                    retained = min(rep_ids, key=lambda t: extract_transcript_version(t))
                else:
                    # No representative; retain the one with the lowest version number.
                    retained = min(mrna_ids, key=lambda t: extract_transcript_version(t))
                # Mark all other transcripts for removal.
                for mrna_id in mrna_ids:
                    if mrna_id != retained:
                        removal_set.add(mrna_id)
                if verbose:
                    duplicate_ids = [mrna_id for mrna_id in mrna_ids if mrna_id != retained]
                    print(f"Duplicate mRNA structure detected in gene {gene_id}:")
                    print(f"  Exon structure: {structure}")
                    print(f"  Retaining transcript: {retained}")
                    print(f"  Removing transcripts: {', '.join(duplicate_ids)}")

        if removal_set:
            # Filter out features associated with mRNAs in the removal set.
            new_features = []
            for f in features:
                # If this is an mRNA (or ncRNA) feature marked for removal, skip it.
                if f['type'] in ['mRNA', 'ncRNA'] and f['id'] in removal_set:
                    continue

                # For child features (exon, CDS, UTR, etc.), check if they are linked to a removed mRNA.
                remove_feature = False
                # Check the parsed 'transcript' field.
                if f.get('transcript') and f['transcript'] in removal_set:
                    remove_feature = True
                # Also check the attributes string for a "Parent=" field.
                elif "Parent=" in f['attributes']:
                    # Example attribute: "ID=...;Parent=TraesCAD_EIv1.0_0468320.2"
                    parent_parts = [part for part in f['attributes'].split(';') if part.startswith("Parent=")]
                    if parent_parts:
                        parent_ids = parent_parts[0].split('=')[1].split(',')
                        for pid in parent_ids:
                            if pid in removal_set:
                                remove_feature = True
                                break
                if remove_feature:
                    continue

                new_features.append(f)
            dictionary[gene_id] = new_features

    return dictionary

def create_new_gff(outfile, edited_dictionary, header, abutting_removals, verbose):
    """
    Write a new GFF3 file with the edited features, ensuring models with abutting exons are removed,
    recalculated exons and their corresponding CDS features are preserved.
    Deduplicate any duplicate CDS entries.
    """
    try:
        if os.path.isfile(outfile):
            raise TypeError("File already exists. Please delete the existing file or provide an alternative outfile name")

        with open(outfile, 'w') as file:
            # Write the GFF3 header
            for line in header:
                file.write(line)

            # Sort genes by seqid and start coordinate
            sorted_genes = sorted(
                edited_dictionary.items(),
                key=lambda x: (x[1][0]['seqid'], min(f['start'] for f in x[1] if f['type'] in ['gene', 'ncRNA_gene']))
            )

            for gene_id, features in sorted_genes:
                # Exclude features associated with transcripts flagged for removal
                filtered_features = [
                    f for f in features if not (
                        f['type'] in ['mRNA', 'ncRNA'] and f['id'] in abutting_removals or
                        f.get('transcript') in abutting_removals
                    )
                ]

                if not filtered_features:  # Skip genes with no valid features remaining
                    continue

                # Separate features into categories
                genes = [f for f in filtered_features if f['type'] in ['gene', 'ncRNA_gene']]
                mrnas = [f for f in filtered_features if f['type'] in ['mRNA', 'ncRNA']]
                child_features = [f for f in filtered_features if f['type'] in ['exon', 'CDS', 'three_prime_UTR', 'five_prime_UTR']]
                recalculated_exons = [f for f in features if f['type'] == 'exon' and 'recalc_exon' in f['attributes']]
                recalculated_cds = [f for f in features if f['type'] == 'CDS' and 'CDS' in f['attributes']]

                # Deduplication mechanism
                written_features = set()

                # Write genes and their child features
                for gene in genes:
                    file.write('\t'.join([
                        gene['seqid'], gene['source'], gene['type'], str(gene['start']), str(gene['end']),
                        gene['score'], gene['strand'], gene['phase'], gene['attributes']
                    ]) + '\n')

                    # Sort mRNAs by start coordinate (they are already grouped by seqid due to parent sort)
                    sorted_mrnas = sorted(mrnas, key=lambda m: m['start'])

                    for mrna in sorted_mrnas:
                        file.write('\t'.join([
                            mrna['seqid'], mrna['source'], mrna['type'], str(mrna['start']), str(mrna['end']),
                            mrna['score'], mrna['strand'], mrna['phase'], mrna['attributes']
                        ]) + '\n')

                        # Collect child features for this mRNA
                        mrna_id = mrna['id']
                        associated_features = [
                            f for f in child_features if f.get('transcript') == mrna_id
                        ] + [
                            f for f in recalculated_exons if f['attributes'].split(';')[1].split('=')[1] == mrna_id
                        ] + [
                            f for f in recalculated_cds if f['attributes'].split(';')[1].split('=')[1] == mrna_id
                        ]

                        # Sort child features by start, end, and feature type
                        sorted_children = sorted(
                            associated_features,
                            key=lambda c: (
                                c['start'],  # Primary sort by start coordinate
                                c['end'],    # Secondary sort by end coordinate
                                c['type']    # Tertiary sort alphabetically by feature type
                            )
                        )

                        # Write the sorted child features
                        for child in sorted_children:
                            # Create a unique identifier for the feature to prevent duplicates
                            feature_id = (child['seqid'], child['start'], child['end'], child['strand'], child['type'], child['attributes'])
                            if feature_id not in written_features:
                                written_features.add(feature_id)
                                file.write('\t'.join([
                                    child['seqid'], child['source'], child['type'], str(child['start']), str(child['end']),
                                    child['score'], child['strand'], child['phase'], child['attributes']
                                ]) + '\n')

                # Write gene separator
                file.write('###\n')

        print(f"New GFF file written to: {outfile}")

    except TypeError as e:
        print(e)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Identify and resolve duplicate fnd abutting eatures in a GFF3 file.")
    parser.add_argument("gff_file", help="Path to the GFF file.")
    parser.add_argument("-o", "--output", default="duplicate_removed.gff3", help="Output file (default: duplicate_removed.gff3).")
    parser.add_argument("-v", "--verbose", action="store_true", default=False, help="Enable verbose output.")
    args = parser.parse_args()

    dictionaries = input_gff_file(args.gff_file)
    HEADER = dictionaries['header']
    GFF3_DICT = dictionaries['gff3_dict']
    EXON = dictionaries['exon']

    # Step 2: Identify duplicate features
    duplicate_genes = get_duplicate_features(GFF3_DICT, args.verbose)

#    if args.verbose:
#        print(f"Duplicated features detected in the following genes: {', '.join(duplicate_genes)}")

    # Step 3: Resolve duplicate features
    final_edited_dictionary = resolve_duplicates_with_cds_based_adjustment(GFF3_DICT, duplicate_genes, args.verbose)

    # Step 4: Identify and resolve abutting exons
    abutting_removals = check_abutting_exons(EXON, args.verbose)
    if args.verbose and abutting_removals:
        print(f"Abutting exons require removal of the following transcripts: {', '.join(abutting_removals)}")

    # Step 5: Remove duplicate mRNAs with identical exon structures.
    final_edited_dictionary = remove_duplicate_mrna_structures(final_edited_dictionary, args.verbose)

    # Step 6: Write the edited dictionary to a new GFF file.
    create_new_gff(args.output, final_edited_dictionary, HEADER, abutting_removals, args.verbose)


