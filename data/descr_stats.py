#!/usr/bin/python
#
# This script generates summary statistics and raw plots for the data note
# associated with the annotations of portrayed emotions in the movie
# Forrest Gump. It is intended to serve as a more detailed description
# of the employed analysis and aggregation procedures than what is possible
# to convey in a manuscript.
#
# In order to reproduce the results, the script needs to be executed in the
# root of the extracted dataset.  Summary statistics are printed with LaTeX
# markup and were directly included into the LaTeX sources of the associated
# Data note publication.
#
# Required arguments:
#   1. path to store the generated inter-observer agreement times series
#   2. path to store the generated figures
#
# The following python packages are required:
#   - NumPy
#   - SciPy
#   - PyMVPA
#   - seaborn
#
# Example:
#   $ python descr_stats.py /tmp /tmp
#
# This source code is (C) by Michael Hanke <michael.hanke@gmail.com> and
# made available under the terms of the Creative Common Attribution-ShareAlike
# 4.0 International (CC BY-SA 4.0) license.
#

import numpy as np
from scipy.stats import spearmanr

# hard code the max duration of the movie stimulus
maxmovietime = 7085.28

#
# Load data
#

def get_shots():
    starts = np.loadtxt('movie_shots.csv')
    segments = np.array((starts,
                         np.concatenate((starts[1:],
                                         (maxmovietime,))))).T
    return segments

def get_scenes():
    starts = np.recfromcsv('movie_scenes.csv',
                           names=('start','title','tod','set'))['start']
    segments = np.array((starts,
                         np.concatenate((starts[1:],
                                         (maxmovietime,))))).T
    return segments

def get_nsecond_segments(n=1):
    max = get_scenes()[-1,1]
    return np.array((np.arange(0,max-n,n), np.arange(n,max,n))).T

def get_av_ratings():
    import glob
    return [np.recfromcsv(f) for f in glob.glob('raw/av*.csv')]

def get_ao_ratings():
    import glob
    return [np.recfromcsv(f) for f in glob.glob('raw/ao*.csv')]

def get_all_ratings():
    return get_av_ratings() + get_ao_ratings()

#
# Stats
#

def get_labeled_fraction(rat, col):
    # what fraction of the annotations carry values in a specific column
    tot = np.sum([len(r) for r in rat])
    lbl = np.sum([len(r) for r in get_labeled_ratings(rat, col)])
    return float(lbl) / tot

def get_agreed_labels(ratings, col, segments, athresh=0.5, nseg_thresh=0):
    # determine values for a particular column that show a minimum
    # inter-observer agreement for a minimum number of time segments
    # anywhere in the movie
    from scipy.ndimage.measurements import label
    labels = \
        np.unique(
            np.concatenate(
                [np.unique(
                    np.concatenate(
                        [d.split() for d in r[col]]))
                    for r in ratings]))
    found = []
    for l in labels:
        match = slice2segments(ratings, {col: l}, segments) > athresh
        nseg = np.sum(match)
        nblobs = label(match)[1]
        if nblobs > nseg_thresh:
            found.append((l, nseg, nblobs))
    return found

def calc_bootstrapped_intercorrelation(ratings, cond1, cond2, segments):
    # split the set of observers into all possible ~halves and
    # compute the time series correlations of inter-oberserver
    # agreement wrt annotation matching particular criteria across
    # both groups
    from mvpa2.misc.support import unique_combinations
    N = len(ratings)
    corr = []
    for combo in unique_combinations(range(N), N / 2):
        half1 = [ratings[i] for i in combo]
        half2 = [ratings[i] for i in xrange(N) if not i in combo]
        c1 = slice2segments(half1, cond1, segments) \
                - slice2segments(half1, cond2, segments)
        c2 = slice2segments(half2, cond1, segments) \
                - slice2segments(half2, cond2, segments)
        corr.append(spearmanr(c1, c2)[0])
    return corr

def get_ci_stats(arr):
    # convert an array of correlation scores into a LaTeX
    # markup with mean and CI
    m = np.mean(arr)
    sem = np.std(arr) / np.sqrt(len(arr))
    if np.isnan(m):
        return 'n/a'
    else:
        if m >= 0.5:
            return '\\textbf{%.3f} $\\pm$%.3f' % (m, 1.96 * sem)
        else:
            return '%.3f $\\pm$%.3f' % (m, 1.96 * sem)

def get_corr_ci(v1, v2):
    # take to time series, compute the correlation, and yield a
    # LaTeX markup of the value plus a CI (via Fisher transform.)
    c = spearmanr(v1, v2)[0]
    # fisher transform
    fc = np.arctanh(c)
    se = 1. / np.sqrt(len(v1) - 3)
    # back to correlation
    ci = np.tanh(fc + 1.96 * se)
    if np.isnan(c):
        return 'n/a'
    else:
        if c >= 0.5:
            return '\\textbf{%.3f} $\\pm$%.3f' % (c, ci - c)
        else:
            return '%.3f $\\pm$%.3f' % (c, ci - c)

def print_stats(rat, rat_label, all_rat):
    # compute various annotation statistics
    print '\\newcommand{\\%sTotalRaters}{%i}' % (rat_label, len(rat))
    athresh = 0.5
    print '\\newcommand{\\%sAggThresh}{%i\\%%}' % (rat_label, athresh * 100)
    segments = get_nsecond_segments()
    print '\\newcommand{\\%sFracWithLabeledChar}{%.1f\\%%}' \
            % (rat_label, get_labeled_fraction(rat, 'character') * 100)
    e = get_agreed_labels(rat, 'character', segments, athresh=-1)
    print '%% %s total character labels' % (rat_label,)
    print '%% %s' % [v[0] for v in e]
    print '\\newcommand{\\%sTotalCharLabels}{%i}' % (rat_label, len(e))
    e = get_agreed_labels(rat, 'character', segments, athresh=athresh, nseg_thresh=5)
    print '%% %s character labels AGG > %.2f' % (rat_label, athresh)
    print '%% %s' % e
    print '\\newcommand{\\%sThreshCharLabels}{%i}' % (rat_label, len(e))
    print '\\newcommand{\\%sFracWithLabeledEmotions}{%.1f\\%%}' \
            % (rat_label, get_labeled_fraction(rat, 'emotion') * 100)
    e = get_agreed_labels(rat, 'emotion', segments, athresh=athresh)
    print '%% %s emotion labels AGG > %.2f' % (rat_label, athresh)
    print '%% %s' % e
    print '\\newcommand{\\%sThreshEmoLabels}{%i}' % (rat_label, len(e))
    print '\\newcommand{\\%sFracWithLabeledOncue}{%.1f\\%%}' \
            % (rat_label, get_labeled_fraction(rat, 'oncue') * 100)
    e = get_agreed_labels(rat, 'oncue', segments, athresh=athresh)
    print '%% %s oncue labels AGG > %.2f' % (rat_label, athresh)
    print '%% %s' % e
    print '\\newcommand{\\%sThreshOncueLabels}{%i}' % (rat_label, len(e))
    print '\\newcommand{\\%sFracWithLabeledOffcue}{%.1f\\%%}' \
            % (rat_label, get_labeled_fraction(rat, 'offcue') * 100)
    e = get_agreed_labels(rat, 'offcue', segments, athresh=athresh)
    print '%% %s offcue labels AGG > %.2f' % (rat_label, athresh)
    print '%% %s' % e
    print '\\newcommand{\\%sThreshOffcueLabels}{%i}' % (rat_label, len(e))

    # per character stats
    for char, clabel in (('*', 'AllChar'),
                         ('FORREST', 'Forrest'),
                         ('JENNY', 'Jenny')):
        print '\\newcommand{\\%sCorrArousalValence%s}{%s}' \
                % (rat_label, clabel,
                   get_corr_ci(get_arousal_modulation(rat, segments, char=char),
                               get_valence_modulation(rat, segments, char=char)))
        print '\\newcommand{\\%sCorrValenceDirection%s}{%s}' \
                % (rat_label, clabel,
                  get_corr_ci(get_valence_modulation(rat, segments, char=char),
                               get_direction_modulation(rat, segments, char=char)))
        print '\\newcommand{\\%sCorrArousalDirection%s}{%s}' \
                % (rat_label, clabel,
                   get_corr_ci(get_arousal_modulation(rat, segments, char=char),
                               get_direction_modulation(rat, segments, char=char)))

        s = get_ci_stats(
                calc_bootstrapped_intercorrelation(
                    rat, {'arousal': 'HIGH', 'character':char},
                         {'arousal':'LOW', 'character': char},
                    segments))
        print '\\newcommand{\\%sInterRaterConsistArousal%s}{%s}' \
                % (rat_label, clabel, s)
        s = get_ci_stats(
                calc_bootstrapped_intercorrelation(
                    rat, {'valence': 'POS', 'character':char},
                         {'valence':'NEG', 'character': char},
                    segments))
        print '\\newcommand{\\%sInterRaterConsistValence%s}{%s}' \
                % (rat_label, clabel, s)
        s = get_ci_stats(
                calc_bootstrapped_intercorrelation(
                    rat, {'direction': 'SELF', 'character':char},
                         {'direction':'OTHER', 'character': char},
                    segments))
        print '\\newcommand{\\%sInterRaterConsistDirection%s}{%s}' \
                % (rat_label, clabel, s)


        for emo in get_unique_emotions(all_rat):
            s = get_ci_stats(
                    calc_bootstrapped_intercorrelation(
                        rat, {'emotion': emo, 'character':char},
                             {'emotion': None, 'character': char},
                        segments))
            print '\\newcommand{\\%sInterRaterConsist%s%s}{%s}' \
                    % (rat_label, emo, clabel, s)
        for cue in get_unique_oncues(all_rat):
            s = get_ci_stats(
                    calc_bootstrapped_intercorrelation(
                        rat, {'oncue': cue, 'character':char},
                             {'oncue': None, 'character': char},
                        segments))
            print '\\newcommand{\\%sInterRaterConsist%s%s}{%s}' \
                    % (rat_label, cue, clabel, s)


#
# Plots
#

def comp_barplot(av, ao, col, ylabel):
    # helper for all bar plots
    import pylab as pl
    from mvpa2.misc.plot import plot_bars
    avd = {e[0]: e[col] for e in av}
    aod = {e[0]: e[col] for e in ao}
    joind = {k: 0 for k in avd.keys() + aod.keys()}
    for k, v in avd.iteritems():
        joind[k] = v
    p = plot_bars([joind[k] for k in sorted(joind.keys())], yerr=None,
              offset=-.15, width=.3)
    p.set_label("audio-visual")
    joind = {k: 0 for k in avd.keys() + aod.keys()}
    for k, v in aod.iteritems():
        joind[k] = v
    p = plot_bars([joind[k] for k in sorted(joind.keys())], yerr=None,
              color='0.2', offset=.15, width=.3,
              labels=[bn(n) for n in sorted(joind.keys())])
    p.set_label("audio-only")
    pl.xlim((-1, len(joind)))
    pl.ylabel(ylabel)
    pl.legend()

def bipolar_tsplot(av, ao, col, colx_labels, ylabels, segments, char='*'):
    # helper for all time series plot for bipolar variables
    ts = slice2segments(ao, {col:colx_labels[0], 'character':char}, segments) \
            - slice2segments(ao, {col:colx_labels[1], 'character':char}, segments)
    pl.plot(segments.T[0], ts, label='audio-only', color='blue', alpha=.5)
    pl.fill_between(segments.T[0], ts, color='blue', alpha=.5)
    ts = slice2segments(av, {col:colx_labels[0], 'character':char}, segments) \
            - slice2segments(av, {col:colx_labels[1], 'character':char}, segments)
    pl.plot(segments.T[0], ts, label='audio-visual', color='green', alpha=.5)
    pl.fill_between(segments.T[0], ts, color='green', alpha=.5)
    pl.ylabel('%s' % (bn(col)))
    pl.xlabel('movie time in seconds')
    pl.xlim((segments[0,0], segments[-1,1]))
    pl.yticks(np.array((-1,0,1)), ylabels)
    pl.ylim((-1,1))
    pl.legend()

def unipolar_tsplot(av, ao, col, colx_label, segments, char='*'):
    # helper for all time series plot for unipolar variables
    pl.fill_between(segments.T[0], slice2segments(ao, {col:colx_label, 'character':char}, segments),
            label='audio-only', zorder=0, color='blue', alpha=.5)
    pl.fill_between(segments.T[0], slice2segments(av, {col:colx_label, 'character':char}, segments),
            label='audio-visual', zorder=1, color='green', alpha=.5)
    pl.ylabel('%s' % (bn(colx_label)))
    pl.xlabel('movie time in seconds')
    pl.xlim((segments[0,0], segments[-1,1]))
    pl.ylim((0,1))
    pl.yticks(np.array((0,1)), ('absent (0)', 'present (1)'))

def mkplot_indicator_ts(avr, aor, character, segments, figpath):
    # demo plot with time series of various properties across the movie
    fig = pl.figure(figsize=(10,8), dpi=300)
    ax=pl.subplot(811)
    bipolar_tsplot(avr, aor,
                   'arousal', ('HIGH', 'LOW'),
                   ('low (-1)', 'neutral (0)', 'high (+1)'),
                   segments,
                   char=character)
    ax=pl.subplot(812)
    bipolar_tsplot(avr, aor,
                   'valence', ('POS', 'NEG'),
                   ('neg (-1)', 'neutral (0)', 'pos (+1)'),
                   segments,
                   char=character)
    ax=pl.subplot(813)
    unipolar_tsplot(avr, aor, 'emotion', 'HAPPINESS', segments, char=character)
    ax=pl.subplot(814)
    unipolar_tsplot(avr, aor, 'emotion', 'LOVE', segments, char=character)
    ax=pl.subplot(815)
    unipolar_tsplot(avr, aor, 'emotion', 'FEAR', segments, char=character)
    ax=pl.subplot(816)
    unipolar_tsplot(avr, aor, 'emotion', 'SADNESS', segments, char=character)
    ax=pl.subplot(817)
    unipolar_tsplot(avr, aor, 'emotion', 'ANGERRAGE', segments, char=character)
    ax=pl.subplot(818)
    unipolar_tsplot(avr, aor, 'oncue', 'VERBAL', segments, char=character)
    fig.autofmt_xdate()
    if character == '*':
        character = 'allchar'
    pl.savefig(opj(figpath, 'indicator_ts_%s.svg' % character.lower()))

def print_combstats(avr, aor):
    # various stats computed across stimulus types
    segments = get_nsecond_segments()
    for char, clabel in (('*', 'AllChar'),
                         ('FORREST', 'Forrest'),
                         ('JENNY', 'Jenny')):
        print '\\newcommand{\\InterModCorrArousal%s}{%s}' \
                % (clabel,
                   get_corr_ci(get_arousal_modulation(avr, segments, char=char),
                              get_arousal_modulation(aor, segments, char=char)))
        print '\\newcommand{\\InterModCorrValence%s}{%s}' \
                % (clabel,
                   get_corr_ci(get_valence_modulation(avr, segments, char=char),
                              get_valence_modulation(aor, segments, char=char)))
        print '\\newcommand{\\InterModCorrDirection%s}{%s}' \
                % (clabel,
                   get_corr_ci(get_direction_modulation(avr, segments, char=char),
                              get_direction_modulation(aor, segments, char=char)))
        for emo in get_unique_emotions(avr + aor):
            print '\\newcommand{\\InterModCorr%s%s}{%s}' \
                    % (emo, clabel,
                       get_corr_ci(_get_modulation(avr, segments, emotion=emo, character=char),
                                   _get_modulation(aor, segments, emotion=emo, character=char)))
        for cue in get_unique_oncues(avr + aor):
            print '\\newcommand{\\InterModCorr%s%s}{%s}' \
                    % (cue, clabel,
                       get_corr_ci(_get_modulation(avr, segments, oncue=cue, character=char),
                                   _get_modulation(aor, segments, oncue=cue, character=char)))
#
# Segmentation
#

def mk_thresh_emotion_episodes(rat, thresh, segments):
    # yield per character list of emotion episodes with a minimum inter-observer
    # agreement wrt any emotion attribute
    chars = get_unique_characters(rat)
    episodes = {}

    def _postprocess(e):
        return {k: np.median(v) for k, v in e.iteritems()}

    for char in chars:
        ep = episodes.get(char, [])
        ind = [get_arousal_modulation(rat, segments, char=char)]
        labels = ['arousal']
        for l, d in (('v_pos', dict(valence='POS')),
                     ('v_neg', dict(valence='NEG')),
                     ('d_self', dict(direction='SELF')),
                     ('d_other', dict(direction='OTHER')),
                     ('e_admiration', dict(emotion='ADMIRATION')),
                     ('e_anger/rage', dict(emotion='ANGER/RAGE')),
                     ('e_contempt',   dict(emotion='CONTEMPT')),
                     ('e_disappointment', dict(emotion='DISAPPOINTMENT')),
                     ('e_fear',       dict(emotion='FEAR')),
                     ('e_fears_confirmed', dict(emotion='FEARS_CONFIRMED')),
                     ('e_gloating',   dict(emotion= 'GLOATING')),
                     ('e_gratification',  dict(emotion='GRATIFICATION')),
                     ('e_gratitude',  dict(emotion='GRATITUDE')),
                     ('e_happiness',  dict(emotion='HAPPINESS')),
                     ('e_happy-for',  dict(emotion='HAPPY-FOR')),
                     ('e_hate', dict(emotion='HATE')),
                     ('e_hope', dict(emotion='HOPE')),
                     ('e_love', dict(emotion='LOVE')),
                     ('e_pity/compassion', dict(emotion= 'PITY/COMPASSION')),
                     ('e_pride',   dict(emotion='PRIDE')),
                     ('e_relief',  dict(emotion='RELIEF')),
                     ('e_remorse', dict(emotion='REMORSE')),
                     ('e_resent',  dict(emotion='RESENTMENT')),
                     ('e_sadness', dict(emotion='SADNESS')),
                     ('e_satisfaction', dict(emotion='SATISFACTION')),
                     ('e_shame',   dict(emotion='SHAME')),
                     ('c_audio',   dict(oncue='AUDIO')),
                     ('c_context',   dict(oncue='CONTEXT')),
                     ('c_face',   dict(oncue='FACE')),
                     ('c_gesture',   dict(oncue='GESTURE')),
                     ('c_narrator',   dict(oncue='NARRATOR')),
                     ('c_verbal',   dict(oncue='VERBAL')),
                    ):
            ind.append(_get_modulation(rat, segments, character=char, **d))
            labels.append(l)
        ind = np.array(ind)
        # where is any above threshold agreement
        flags = np.abs(ind) >= thresh
        staging = None
        last_ind = np.array([False] * len(ind))
        # for each segment
        for i, f in enumerate(flags.T):
            #print i, f,
            if not np.sum(f):
                if staging:
                    ep.append(_postprocess(staging))
                    staging = None
                    #print 'commit',
                last_ind = f
                #print 'skip'
                continue
            # continuing episode?
            if np.all(f == last_ind):
                # end of annotation is end of current segment
                staging['end'] = segments[i, 1]
                for nl, l in enumerate(labels):
                    staging[l].append(ind[nl, i])
                #print 'extend'
            else:
                # new episode
                if staging:
                    #print 'commit',
                    ep.append(_postprocess(staging))
                #print 'new'
                staging = dict(start=segments[i, 0],
                               end=segments[i, 1])
                last_ind = f
                for nl, l in enumerate(labels):
                    staging[l] = [ind[nl, i]]

        episodes[char] = ep
    return episodes, labels

def emo2advene(data, labels, thresh=0.5):
    # format output of `mk_thresh_emotion_episodes()` into a format that is
    # importable by Advene, while merging all episodes of all characters
    # into a single file
    episodes = []
    s = ''
    for char, ep in data.iteritems():
        for e in ep:
            e['character'] = char
            episodes.append(e)
    episodes = sorted(episodes, cmp=lambda x,y: cmp(x['start'], y['start']))
    for e in episodes:
        tags = []
        if e['arousal'] > thresh:
            tags.append('ha')
        elif e['arousal'] < (-1 * thresh):
            tags.append('la')
        for l in labels:
            if l == 'arousal':
                continue
            if e[l] > thresh:
                tags.append(l[2:])
        e['tags'] = ','.join(tags)
        s += "%(start).1f\t%(end).1f\tchar=%(character)s tags=%(tags)s arousal=%(arousal).2f val_pos=%(v_pos).2f val_neg=%(v_neg).2f\n" % e
    return s

#
# Helpers
#

def bn(n):
    # beautify names
    n = n.lower()
    if n == 'forrestvo':
        return 'Forrest (VO)'
    elif n == 'mrsgump':
        return 'Mrs. Gump'
    elif n == 'disappointment':
        return 'Disappoint.'
    elif n == 'angerrage':
        return 'Anger/Rage'
    return n.capitalize()


def get_unique_characters(rat):
    return np.unique(
            np.concatenate(
                [np.unique([a['character'] for a in an])
                    for an in rat]))

def get_unique_emotions(rat):
    return [e for e in np.unique(
            np.concatenate(
                [np.unique(
                    np.concatenate([a['emotion'].split() for a in an]))
                    for an in rat])) if not '?' in e]

def get_unique_oncues(rat):
    return [e for e in np.unique(
            np.concatenate(
                [np.unique(
                    np.concatenate([a['oncue'].split() for a in an]))
                    for an in rat])) if not '?' in e]


def slice2segments(ratings, cond, segments):
    # compute a time series of inter-observer agreement wrt a particular
    # emotion property (or combinations thereof) 
    # annotations given with start and stop time, are converted into a
    # timeseries with data point locations given by the sequence of
    # `segments`. Segments intersecting with a given annotation from an
    # individual observer are set to one, the rest to zero. The mean
    # across observers for any segment is returned
    slicer = np.zeros(len(segments))
    for rat in ratings:
        rslicer = np.zeros(len(segments))
        for e in rat:
            use = True
            for k, v in cond.iteritems():
                if v == '*':
                    continue
                if k in ('oncue', 'offcue', 'emotion'):
                    if not v in e[k].split():
                        use = False
                else:
                    if not v == e[k]:
                        use = False
            if not use:
                continue
            select = np.logical_and(segments.T[1] > e['start'],
                                    segments.T[0] < e['end'])
            rslicer[select] += 1
        slicer += rslicer > 0
    slicer = slicer.astype(float) / len(ratings)
    return slicer

def get_labeled_ratings(rat, col):
    return [r[r[col] != ''] for r in rat]

def get_timeseries(rat, urat, segments, char='*'):
    # yield time series representations of all relevant emotion attributes
    # from raw annotations
    vars = [get_arousal_modulation(rat, segments, char=char),
            get_valence_modulation(rat, segments, char=char),
            get_direction_modulation(rat, segments, char=char)]
    labels = ['arousal', 'valence', 'direction']
    for emo in get_unique_emotions(urat):
        vars.append(_get_modulation(rat, segments, emotion=emo, character=char))
        labels.append(emo.lower())
    for oc in get_unique_oncues(urat):
        vars.append(_get_modulation(rat, segments, oncue=oc, character=char))
        labels.append(oc.lower())
    return np.array(vars).T, labels

def _get_modulation(ratings, segments, **kwargs):
    return slice2segments(ratings, kwargs, segments)

def get_arousal_modulation(ratings, segments, char='*'):
    ts = _get_modulation(ratings, segments, character=char, arousal='HIGH') \
         - _get_modulation(ratings, segments, character=char, arousal='LOW')
    return ts

def get_valence_modulation(ratings, segments, char='*'):
    ts = _get_modulation(ratings, segments, character=char, valence='POS') \
         - _get_modulation(ratings, segments, character=char, valence='NEG')
    return ts

def get_direction_modulation(ratings, segments, char='*'):
    ts = _get_modulation(ratings, segments, character=char, direction='SELF') \
         - _get_modulation(ratings, segments, character=char, direction='OTHER')
    return ts

if __name__ == '__main__':
    # main function: compute stats, generate derived data, make figures
    import sys
    import pylab as pl
    import os
    from os.path import join as opj

    if len(sys.argv) < 3:
        print "insufficient number of arguments"
        print "usage: descr_stats.py <stats export path> <figure export path>"
        sys.exit(1)

    statspath=sys.argv[1]
    figpath=sys.argv[2]

    for p in (opj(statspath, 'timeseries'),
              opj(statspath, 'segmentation')):
        if not os.path.exists(p):
            os.makedirs(p)

    second_segments = get_nsecond_segments()

    avr = get_av_ratings()
    aor = get_ao_ratings()

    open(opj(statspath, 'segmentation', 'emotions_av_1s_thr50.tsv'), 'w').write(
            emo2advene(
                *mk_thresh_emotion_episodes(avr, .5, get_nsecond_segments(1)),
                thresh=.5))
    open(opj(statspath, 'segmentation', 'emotions_ao_1s_thr50.tsv'), 'w').write(
            emo2advene(
                *mk_thresh_emotion_episodes(aor, .5, get_nsecond_segments(1)),
                thresh=.5))
    open(opj(statspath, 'segmentation', 'emotions_av_shots_thr50.tsv'), 'w').write(
            emo2advene(
                *mk_thresh_emotion_episodes(avr, .5, get_shots()),
                thresh=.5))
    open(opj(statspath, 'segmentation', 'emotions_ao_shots_thr50.tsv'), 'w').write(
            emo2advene(
                *mk_thresh_emotion_episodes(aor, .5, get_shots()),
                thresh=.5))
    # export inter-observer agreement timeseries
    # for stim type
    for rat, ratlabel in ((avr, 'av'), (aor, 'ao')):
        # for segment type
        for seg, seglabel in ((get_nsecond_segments(1), '1s'),
                              (get_nsecond_segments(2), '2s'),
                              (get_shots(), 'shots')):
            # for all characters
            for char in ('*',) + tuple(get_unique_characters(rat)):
                arr, l = get_timeseries(rat, rat, seg, char=char)
                if char == '*':
                    char = 'allchar'
                np.savetxt(
                    opj(statspath, 'timeseries',
                        'ioats_%s_%s_%s.csv' % (seglabel, ratlabel, char.lower())),
                           arr, fmt='%.2f', delimiter=',', comments = '',
                           header=','.join(l))
    # export summary stats for the data paper
    print_stats(avr, 'AV', avr+aor)
    print_stats(aor, 'AO', avr+aor)
    print_combstats(avr, aor)

    # IOA time series demo plot
    mkplot_indicator_ts(avr, aor, '*', get_nsecond_segments(n=10), figpath)

    # episodes per character plot
    fig = pl.figure(figsize=(10,6))
    pl.subplot(121)
    comp_barplot(
        get_agreed_labels(avr, 'character', second_segments, athresh=.5, nseg_thresh=5),
        get_agreed_labels(aor, 'character', second_segments, athresh=.5, nseg_thresh=5),
        2,
        'Total number of emotion episodes')
    ax=pl.subplot(122)
    comp_barplot(
        get_agreed_labels(avr, 'character', second_segments, athresh=.5, nseg_thresh=5),
        get_agreed_labels(aor, 'character', second_segments, athresh=.5, nseg_thresh=5),
        1,
        'Total time of portrayed emotions (in sec)')
    ax.yaxis.tick_right()
    ax.yaxis.set_label_position("right")
    fig.autofmt_xdate()
    pl.savefig(opj(figpath, 'character_episodes.svg'))

    # episodes per emotion catgeory plot
    fig = pl.figure(figsize=(10,6))
    pl.subplot(121)
    comp_barplot(
        get_agreed_labels(avr, 'emotion', second_segments, athresh=.5, nseg_thresh=5),
        get_agreed_labels(aor, 'emotion', second_segments, athresh=.5, nseg_thresh=5),
        2,
        'Total number of emotion episodes')
    ax=pl.subplot(122)
    comp_barplot(
        get_agreed_labels(avr, 'emotion', second_segments, athresh=.5, nseg_thresh=5),
        get_agreed_labels(aor, 'emotion', second_segments, athresh=.5, nseg_thresh=5),
        1,
        'Total number of emotion episodes')
    ax.yaxis.tick_right()
    ax.yaxis.set_label_position("right")
    fig.autofmt_xdate()
    pl.savefig(opj(figpath, 'labeledemotion_episodes.svg'))

    # episodes per onset cue type plot
    fig = pl.figure(figsize=(5,5))
    comp_barplot(
        get_agreed_labels(avr, 'oncue', second_segments, athresh=.5, nseg_thresh=5),
        get_agreed_labels(aor, 'oncue', second_segments, athresh=.5, nseg_thresh=5),
        2,
        'Total number of emotion episodes')
    fig.autofmt_xdate()
    pl.savefig(opj(figpath, 'labeledoncue_episodes.svg'))

    # intra-stimulus IOA time series correlation plot
    import seaborn as sns
    fig = pl.figure()
    cmap_range=(-.6,.6)
    pl.subplot(121)
    m, l = get_timeseries(avr, avr+aor, second_segments)
    l = [bn(i) for i in l]
    sns.corrplot(m, annot=False, method='spearman', names=l, diag_names=False, cmap_range=cmap_range)
    pl.subplot(122)
    m, l = get_timeseries(aor, avr+aor, second_segments)
    l = [bn(i) for i in l]
    sns.corrplot(m, annot=False, method='spearman', names=l, diag_names=False, cmap_range=cmap_range)
    #fig.autofmt_xdate()
    pl.savefig(opj(figpath, 'intercorr_indicators.svg'))
    #pl.show()
