# -*- coding: utf-8 -*-
r"""
Graphs from the World Map

The methods defined here appear in :mod:`sage.graphs.graph_generators`.
"""

# ****************************************************************************
#
#           Copyright (C) 2006 Robert L. Miller <rlmillster@gmail.com>
#                              and Emily A. Kirkman
#           Copyright (C) 2009 Michael C. Yurko <myurko@gmail.com>
#
# Distributed  under  the  terms  of  the  GNU  General  Public  License (GPL)
#                         http://www.gnu.org/licenses/
# ****************************************************************************

# import from Sage library
from sage.graphs.graph import Graph

def AfricaMap(continental=False, year=2018):
    """
    Return African states as a graph of common border.

    "African state" here is defined as an independent state having the capital
    city in Africa. The graph has an edge between those countries that have
    common *land* border.

    INPUT:

    - ``continental`` -- boolean (default: ``False``); whether to only return
      states in the continental Africa or all African states

    - ``year`` -- integer (default: ``2018``); reserved for future use

    EXAMPLES::

        sage: Africa = graphs.AfricaMap(); Africa
        Africa Map: Graph on 54 vertices
        sage: sorted(Africa.neighbors('Libya'))
        ['Algeria', 'Chad', 'Egypt', 'Niger', 'Sudan', 'Tunisia']

        sage: cont_Africa = graphs.AfricaMap(continental=True)
        sage: cont_Africa.order()
        48
        sage: 'Madagaskar' in cont_Africa
        False

    TESTS::

        sage: Africa.plot()  # long time
        Graphics object consisting of 159 graphics primitives
    """
    if year != 2018:
        raise ValueError("currently only year 2018 is implemented")

    common_border = {
     'Algeria': ['Libya', 'Mali', 'Mauritania', 'Morocco', 'Niger', 'Tunisia'],
     'Angola': ['Namibia', 'Zambia'],
     'Benin': ['Burkina Faso', 'Niger', 'Nigeria', 'Togo'],
     'Botswana': ['Namibia', 'South Africa', 'Zimbabwe'],
     'Burkina Faso': ['Ghana', 'Ivory Coast', 'Mali', 'Niger', 'Togo'],
     'Cameroon': ['Central Africa', 'Chad', 'Equatorial Guinea', 'Gabon',
                      'Nigeria'],
     'Central Africa': ['Chad', 'South Sudan', 'Sudan'],
     'Chad': ['Libya', 'Niger', 'Nigeria', 'Sudan'],
     'Republic of the Congo': ['Gabon', 'Cameroon', 'Central Africa', 'Angola',
                               'Democratic Republic of the Congo'],
     'Democratic Republic of the Congo': ['Zambia', 'South Sudan', 'Tanzania',
                                          'Burundi', 'Rwanda', 'Uganda',
                                          'Central Africa', 'Angola'],
     'Djibouti': ['Eritrea', 'Ethiopia', 'Somalia'],
     'Ethiopia': ['Eritrea', 'Kenya', 'Somalia', 'South Sudan', 'Sudan'],
     'Gabon': ['Equatorial Guinea'],
     'Ghana': ['Ivory Coast', 'Togo'],
     'Guinea': ['Guinea-Bissau', 'Ivory Coast', 'Liberia', 'Sierra Leone'],
     'Kenya': ['Somalia', 'South Sudan', 'Tanzania', 'Uganda'],
     'Liberia': ['Ivory Coast', 'Sierra Leone'],
     'Libya': ['Egypt', 'Niger', 'Sudan', 'Tunisia'],
     'Mali': ['Guinea', 'Ivory Coast', 'Mauritania', 'Niger', 'Senegal'],
     'Mozambique': ['Malawi', 'South Africa', 'Swaziland', 'Zimbabwe'],
     'Niger': ['Nigeria'],
     'Rwanda': ['Burundi', 'Tanzania', 'Uganda'],
     'Senegal': ['Guinea', 'Guinea-Bissau', 'Mauritania', 'Gambia'],
     'South Africa': ['Lesotho', 'Namibia', 'Swaziland', 'Zimbabwe'],
     'South Sudan': ['Uganda', 'Sudan', 'Democratic Republic of the Congo'],
     'Sudan': ['Egypt', 'Eritrea'],
     'Tanzania': ['Burundi', 'Malawi', 'Mozambique', 'Uganda', 'Zambia'],
     'Zambia': ['Malawi', 'Mozambique', 'Namibia', 'Zimbabwe']
     }

    no_land_border = ['Cape Verde', 'Seychelles', 'Mauritius',
                      'São Tomé and Príncipe', 'Madagascar', 'Comoros']

    G = Graph(common_border, format='dict_of_lists')

    if continental:
        G = G.subgraph(G.connected_component_containing_vertex('Central Africa'))
        G.name(new="Continental Africa Map")
    else:
        G.add_vertices(no_land_border)
        G.name(new="Africa Map")

    return G

def EuropeMap(continental=False, year=2018):
    """
    Return European states as a graph of common border.

    "European state" here is defined as an independent state having the capital
    city in Europe. The graph has an edge between those countries that have
    common *land* border.

    INPUT:

    - ``continental`` -- boolean (default: ``False``); whether to only return
      states in the continental Europe or all European states

    - ``year`` -- integer (default: ``2018``); reserved for future use

    EXAMPLES::

        sage: Europe = graphs.EuropeMap(); Europe
        Europe Map: Graph on 44 vertices
        sage: Europe.neighbors('Ireland')
        ['United Kingdom']

        sage: cont_Europe = graphs.EuropeMap(continental=True)
        sage: cont_Europe.order()
        40
        sage: 'Iceland' in cont_Europe
        False
    """
    if year != 2018:
        raise ValueError("currently only year 2018 is implemented")

    common_border = {
     'Austria': ['Czech Republic', 'Germany', 'Liechtenstein', 'Slovenia',
                     'Switzerland'],
     'Belarus': ['Latvia', 'Lithuania', 'Poland', 'Russia', 'Ukraine'],
     'Belgium': ['France', 'Germany', 'Luxembourg', 'Netherlands'],
     'Croatia': ['Bosnia and Herzegovina', 'Hungary', 'Montenegro', 'Serbia',
                     'Slovenia'],
     'France': ['Andorra', 'Germany', 'Italy', 'Luxembourg', 'Monaco',
                    'Switzerland'],
     'Germany': ['Czech Republic', 'Denmark', 'Luxembourg', 'Netherlands',
                     'Switzerland'],
     'Greece': ['Albania', 'Bulgaria', 'Macedonia'],
     'Hungary': ['Austria', 'Romania', 'Serbia', 'Slovakia', 'Slovenia',
                     'Ukraine'],
     'Ireland': ['United Kingdom'],
     'Italy': ['Austria', 'San Marino', 'Slovenia', 'Switzerland',
                   'Vatican City'],
     'Latvia': ['Estonia', 'Lithuania', 'Russia'],
     'Macedonia': ['Albania', 'Bulgaria', 'Serbia'],
     'Montenegro': ['Albania', 'Bosnia and Herzegovina', 'Serbia'],
     'Norway': ['Finland', 'Russia', 'Sweden'],
     'Poland': ['Czech Republic', 'Germany', 'Lithuania', 'Russia', 'Slovakia',
                    'Ukraine'],
     'Romania': ['Bulgaria', 'Moldova', 'Serbia', 'Ukraine'],
     'Russia': ['Estonia', 'Finland', 'Lithuania', 'Ukraine'],
     'Serbia': ['Bosnia and Herzegovina', 'Bulgaria'],
     'Slovakia': ['Austria', 'Czech Republic', 'Ukraine'],
     'Spain': ['Andorra', 'France', 'Portugal'],
     'Sweden': ['Finland'],
     'Switzerland': ['Liechtenstein'],
     'Ukraine': ['Moldova']
    }
    no_land_border = ['Iceland', 'Malta']

    G = Graph(common_border, format='dict_of_lists')

    if continental:
        G = G.subgraph(G.connected_component_containing_vertex('Austria'))
        G.name(new="Continental Europe Map")
    else:
        G.add_vertices(no_land_border)
        G.name(new="Europe Map")

    return G

def USAMap(continental=False):
    """
    Return states of USA as a graph of common border.

    The graph has an edge between those states that have common *land* border
    line or point. Hence for example Colorado and Arizona are marked as
    neighbors, but Michigan and Minnesota are not.

    INPUT:

    - ``continental`` -- boolean (default: ``False``); whether to exclude Alaska
      and Hawaii

    EXAMPLES:

    How many states are neighbor's neighbor for Pennsylvania::

        sage: USA = graphs.USAMap()
        sage: distance = USA.shortest_path_lengths('Pennsylvania')
        sage: len([n2 for n2, d in distance.items() if d == 2])
        7

    Diameter for continental USA::

        sage: USAcont = graphs.USAMap(continental=True)
        sage: USAcont.diameter()
        11
    """
    states = {
    "Alabama": ["Florida", "Georgia", "Mississippi", "Tennessee"],
    "Arizona": ["California", "Colorado", "Nevada", "New Mexico", "Utah"],
    "Arkansas": ["Louisiana", "Mississippi", "Missouri", "Oklahoma",
                     "Tennessee", "Texas"],
    "California": ["Arizona", "Nevada", "Oregon"],
    "Colorado": ["Arizona", "Kansas", "Nebraska", "New Mexico", "Oklahoma",
                     "Utah", "Wyoming"],
    "Connecticut": ["Massachusetts", "New York", "Rhode Island"],
    "Delaware": ["Maryland", "New Jersey", "Pennsylvania"],
    "Florida": ["Alabama", "Georgia"],
    "Georgia": ["Alabama", "Florida", "North Carolina", "South Carolina",
                    "Tennessee"],
    "Idaho": ["Montana", "Nevada", "Oregon", "Utah", "Washington", "Wyoming"],
    "Illinois": ["Indiana", "Iowa", "Michigan", "Kentucky", "Missouri",
                     "Wisconsin"],
    "Indiana": ["Illinois", "Kentucky", "Michigan", "Ohio"],
    "Iowa": ["Illinois", "Minnesota", "Missouri", "Nebraska", "South Dakota",
                 "Wisconsin"],
    "Kansas": ["Colorado", "Missouri", "Nebraska", "Oklahoma"],
    "Kentucky": ["Illinois", "Indiana", "Missouri", "Ohio", "Tennessee",
                     "Virginia", "West Virginia"],
    "Louisiana": ["Arkansas", "Mississippi", "Texas"],
    "Maine": ["New Hampshire"],
    "Maryland": ["Delaware", "Pennsylvania", "Virginia", "West Virginia"],
    "Massachusetts": ["Connecticut", "New Hampshire", "New York",
                          "Rhode Island", "Vermont"],
    "Michigan": ["Illinois", "Indiana", "Ohio", "Wisconsin"],
    "Minnesota": ["Iowa", "North Dakota", "South Dakota", "Wisconsin"],
    "Mississippi": ["Alabama", "Arkansas", "Louisiana", "Tennessee"],
    "Missouri": ["Arkansas", "Illinois", "Iowa", "Kansas", "Kentucky",
                    "Nebraska", "Oklahoma", "Tennessee"],
    "Montana": ["Idaho", "North Dakota", "South Dakota", "Wyoming"],
    "Nebraska": ["Colorado", "Iowa", "Kansas", "Missouri", "South Dakota",
                     "Wyoming"],
    "Nevada": ["Arizona", "California", "Idaho", "Oregon", "Utah"],
    "New Hampshire": ["Maine", "Massachusetts", "Vermont"],
    "New Jersey": ["Delaware", "New York", "Pennsylvania"],
    "New Mexico": ["Arizona", "Colorado", "Oklahoma", "Texas", "Utah"],
    "New York": ["Connecticut", "Massachusetts", "New Jersey",
                     "Pennsylvania", "Vermont"],
    "North Carolina": ["Georgia", "South Carolina", "Tennessee", "Virginia"],
    "North Dakota": ["Minnesota", "Montana", "South Dakota"],
    "Ohio": ["Indiana", "Kentucky", "Michigan", "Pennsylvania",
                 "West Virginia"],
    "Oklahoma": ["Arkansas", "Colorado", "Kansas", "Missouri",
                     "New Mexico", "Texas"],
    "Oregon": ["California", "Idaho", "Nevada", "Washington"],
    "Pennsylvania": ["Delaware", "Maryland", "New Jersey", "New York",
                         "Ohio", "West Virginia"],
    "Rhode Island": ["Connecticut", "Massachusetts"],
    "South Carolina": ["Georgia", "North Carolina"],
    "South Dakota": ["Iowa", "Minnesota", "Montana", "Nebraska",
                         "North Dakota", "Wyoming"],
    "Tennessee": ["Alabama", "Arkansas", "Georgia", "Kentucky", "Mississippi",
                      "Missouri", "North Carolina", "Virginia"],
    "Texas": ["Arkansas", "Louisiana", "New Mexico", "Oklahoma"],
    "Utah": ["Arizona", "Colorado", "Idaho", "Nevada", "New Mexico", "Wyoming"],
    "Vermont": ["Massachusetts", "New Hampshire", "New York"],
    "Virginia": ["Kentucky", "Maryland", "North Carolina", "Tennessee",
                     "West Virginia"],
    "Washington": ["Idaho", "Oregon"],
    "West Virginia": ["Kentucky", "Maryland", "Ohio", "Pennsylvania",
                          "Virginia"],
    "Wisconsin": ["Illinois", "Iowa", "Michigan", "Minnesota"],
    "Wyoming": ["Colorado", "Idaho", "Montana", "Nebraska", "South Dakota",
                    "Utah"]
    }
    if continental:
        name = "Continental USA Map"
    else:
        states['Alaska'] = []
        states['Hawaii'] = []
        name = "USA Map"

    return Graph(states, format='dict_of_lists', name=name)

def WorldMap():
    """
    Return the Graph of all the countries, in which two countries are adjacent
    in the graph if they have a common boundary.

    This graph has been built from the data available
    in The CIA World Factbook [CIA]_ (2009-08-21).

    The returned graph ``G`` has a member ``G.gps_coordinates`` equal to a
    dictionary containing the GPS coordinates of each country's capital city.

    EXAMPLES::

        sage: g = graphs.WorldMap()
        sage: g.has_edge("France", "Italy")
        True
        sage: g.gps_coordinates["Bolivia"]
        [[17, 'S'], [65, 'W']]
        sage: sorted(g.connected_component_containing_vertex('Ireland'))
        ['Ireland', 'United Kingdom']

    TESTS::

    :trac:`24488`::

        sage: 'Iceland' in graphs.WorldMap()
        True
    """
    edges = [
        ('Afghanistan', 'China'), ('Afghanistan', 'Iran'),
        ('Afghanistan', 'Uzbekistan'), ('Albania', 'Greece'),
        ('Albania', 'Kosovo'), ('Albania', 'Macedonia'),
        ('Albania', 'Montenegro'), ('Algeria', 'Morocco'),
        ('Algeria', 'Tunisia'), ('Andorra', 'Spain'),
        ('Angola', 'Democratic Republic of the Congo'), ('Angola', 'Namibia'),
        ('Angola', 'Zambia'), ('Argentina', 'Bolivia'), ('Argentina', 'Brazil'),
        ('Argentina', 'Chile'), ('Argentina', 'Paraguay'),
        ('Argentina', 'Uruguay'), ('Armenia', 'Georgia'), ('Armenia', 'Iran'),
        ('Austria', 'Germany'), ('Azerbaijan', 'Armenia'),
        ('Azerbaijan', 'Georgia'), ('Azerbaijan', 'Iran'),
        ('Azerbaijan', 'Russia'), ('Azerbaijan', 'Turkey'),
        ('Bangladesh', 'Burma'), ('Belgium', 'Germany'),
        ('Belgium', 'Netherlands'), ('Belize', 'Mexico'),
        ('Benin', 'Burkina Faso'), ('Benin', 'Niger'), ('Benin', 'Nigeria'),
        ('Benin', 'Togo'), ('Bolivia', 'Brazil'), ('Bolivia', 'Chile'),
        ('Bolivia', 'Paraguay'), ('Bolivia', 'Peru'),
        ('Bosnia and Herzegovina', 'Croatia'),
        ('Bosnia and Herzegovina', 'Montenegro'),
        ('Bosnia and Herzegovina', 'Serbia'), ('Brazil', 'Colombia'),
        ('Brazil', 'Guyana'), ('Brazil', 'Suriname'), ('Brazil', 'Venezuela'),
        ('Bulgaria', 'Greece'), ('Bulgaria', 'Macedonia'),
        ('Bulgaria', 'Romania'), ('Bulgaria', 'Serbia'),
        ('Burkina Faso', 'Mali'), ('Burkina Faso', 'Niger'),
        ('Burkina Faso', 'Togo'),
        ('Burundi', 'Democratic Republic of the Congo'), ('Cambodia', 'Laos'),
        ('Cambodia', 'Thailand'), ('Cambodia', 'Vietnam'),
        ('Cameroon', 'Central African Republic'), ('Cameroon', 'Chad'),
        ('Cameroon', 'Equatorial Guinea'), ('Cameroon', 'Nigeria'),
        ('Cameroon', 'Republic of the Congo'), ('Canada', 'United States'),
        ('Central African Republic', 'Chad'),
        ('Central African Republic', 'Democratic Republic of the Congo'),
        ('Central African Republic', 'Sudan'), ('Chad', 'Niger'),
        ('Chad', 'Nigeria'), ('Chad', 'Sudan'), ('China', 'Bhutan'),
        ('China', 'Burma'), ('China', 'Hong Kong'), ('China', 'Kazakhstan'),
        ('China', 'Kyrgyzstan'), ('China', 'Mongolia'), ('China', 'Nepal'),
        ('China', 'North Korea'), ('China', 'Russia'), ('China', 'Vietnam'),
        ('Colombia', 'Venezuela'), ('Costa Rica', 'Nicaragua'),
        ("Cote d'Ivoire", 'Burkina Faso'), ("Cote d'Ivoire", 'Guinea'),
        ("Cote d'Ivoire", 'Mali'), ('Cyprus', 'Akrotiri'),
        ('Cyprus', 'Dhekelia'), ('Czech Republic', 'Austria'),
        ('Czech Republic', 'Germany'), ('Czech Republic', 'Poland'),
        ('Democratic Republic of the Congo', 'Zambia'), ('Denmark', 'Germany'),
        ('Djibouti', 'Eritrea'), ('Dominican Republic', 'Haiti'),
        ('Ecuador', 'Colombia'), ('El Salvador', 'Honduras'),
        ('Ethiopia', 'Djibouti'), ('Ethiopia', 'Eritrea'),
        ('Ethiopia', 'Kenya'), ('Ethiopia', 'Somalia'), ('Ethiopia', 'Sudan'),
        ('Finland', 'Russia'), ('Finland', 'Sweden'), ('France', 'Andorra'),
        ('France', 'Belgium'), ('France', 'Brazil'), ('France', 'Germany'),
        ('France', 'Italy'), ('France', 'Luxembourg'), ('France', 'Spain'),
        ('France', 'Suriname'), ('France', 'Switzerland'),
        ('Gabon', 'Cameroon'), ('Gabon', 'Equatorial Guinea'),
        ('Gabon', 'Republic of the Congo'), ('Gaza Strip', 'Egypt'),
        ('Gaza Strip', 'Israel'), ('Ghana', 'Burkina Faso'),
        ('Ghana', "Cote d'Ivoire"), ('Ghana', 'Togo'), ('Gibraltar', 'Spain'),
        ('Guatemala', 'Belize'), ('Guatemala', 'El Salvador'),
        ('Guatemala', 'Honduras'), ('Guatemala', 'Mexico'),
        ('Guinea', 'Sierra Leone'), ('Guinea-Bissau', 'Guinea'),
        ('Guinea-Bissau', 'Senegal'), ('Honduras', 'Nicaragua'),
        ('Hungary', 'Austria'), ('Hungary', 'Croatia'), ('Hungary', 'Serbia'),
        ('India', 'Bangladesh'), ('India', 'Bhutan'), ('India', 'Burma'),
        ('India', 'China'), ('India', 'Nepal'),
        ('Indonesia', 'Papua New Guinea'), ('Iran', 'Iraq'),
        ('Ireland', 'United Kingdom'), ('Israel', 'Egypt'),
        ('Italy', 'Austria'), ('Jordan', 'Iraq'), ('Jordan', 'Israel'),
        ('Jordan', 'Syria'), ('Jordan', 'West Bank'),
        ('Kazakhstan', 'Kyrgyzstan'), ('Kenya', 'Somalia'), ('Kenya', 'Sudan'),
        ('Kenya', 'Uganda'), ('Kosovo', 'Macedonia'), ('Kosovo', 'Serbia'),
        ('Kuwait', 'Iraq'), ('Laos', 'Burma'), ('Laos', 'China'),
        ('Laos', 'Thailand'), ('Laos', 'Vietnam'), ('Latvia', 'Belarus'),
        ('Latvia', 'Estonia'), ('Lebanon', 'Israel'),
        ('Lesotho', 'South Africa'), ('Liberia', "Cote d'Ivoire"),
        ('Liberia', 'Guinea'), ('Liberia', 'Sierra Leone'),
        ('Libya', 'Algeria'), ('Libya', 'Chad'), ('Libya', 'Egypt'),
        ('Libya', 'Niger'), ('Libya', 'Sudan'), ('Libya', 'Tunisia'),
        ('Liechtenstein', 'Austria'), ('Liechtenstein', 'Switzerland'),
        ('Lithuania', 'Belarus'), ('Lithuania', 'Latvia'),
        ('Lithuania', 'Poland'), ('Lithuania', 'Russia'),
        ('Luxembourg', 'Belgium'), ('Luxembourg', 'Germany'),
        ('Macau', 'China'), ('Macedonia', 'Greece'), ('Macedonia', 'Serbia'),
        ('Malaysia', 'Brunei'), ('Malaysia', 'Indonesia'),
        ('Malaysia', 'Thailand'), ('Mali', 'Algeria'), ('Mali', 'Guinea'),
        ('Mali', 'Niger'), ('Mali', 'Senegal'), ('Mauritania', 'Algeria'),
        ('Mauritania', 'Mali'), ('Mauritania', 'Senegal'),
        ('Mauritania', 'Western Sahara'), ('Monaco', 'France'),
        ('Montenegro', 'Croatia'), ('Montenegro', 'Kosovo'),
        ('Montenegro', 'Serbia'), ('Morocco', 'Spain'),
        ('Mozambique', 'Malawi'), ('Mozambique', 'Zambia'),
        ('Mozambique', 'Zimbabwe'), ('Namibia', 'Botswana'),
        ('Namibia', 'Zambia'), ('Netherlands', 'Germany'), ('Niger', 'Algeria'),
        ('Niger', 'Nigeria'), ('Norway', 'Finland'), ('Norway', 'Russia'),
        ('Norway', 'Sweden'), ('Oman', 'United Arab Emirates'),
        ('Oman', 'Yemen'), ('Pakistan', 'Afghanistan'), ('Pakistan', 'China'),
        ('Pakistan', 'India'), ('Pakistan', 'Iran'), ('Panama', 'Colombia'),
        ('Panama', 'Costa Rica'), ('Paraguay', 'Brazil'), ('Peru', 'Brazil'),
        ('Peru', 'Chile'), ('Peru', 'Colombia'), ('Peru', 'Ecuador'),
        ('Poland', 'Belarus'), ('Poland', 'Germany'), ('Portugal', 'Spain'),
        ('Republic of the Congo', 'Angola'),
        ('Republic of the Congo', 'Central African Republic'),
        ('Republic of the Congo', 'Democratic Republic of the Congo'),
        ('Romania', 'Hungary'), ('Romania', 'Moldova'), ('Romania', 'Serbia'),
        ('Russia', 'Belarus'), ('Russia', 'Estonia'), ('Russia', 'Georgia'),
        ('Russia', 'Kazakhstan'), ('Russia', 'Latvia'), ('Russia', 'Mongolia'),
        ('Russia', 'North Korea'), ('Russia', 'Poland'), ('Rwanda', 'Burundi'),
        ('Rwanda', 'Democratic Republic of the Congo'), ('Rwanda', 'Uganda'),
        ('Saint Martin', 'Netherlands Antilles'), ('San Marino', 'Italy'),
        ('Saudi Arabia', 'Iraq'), ('Saudi Arabia', 'Jordan'),
        ('Saudi Arabia', 'Kuwait'), ('Saudi Arabia', 'Oman'),
        ('Saudi Arabia', 'Qatar'), ('Saudi Arabia', 'United Arab Emirates'),
        ('Saudi Arabia', 'Yemen'), ('Senegal', 'Guinea'), ('Serbia', 'Croatia'),
        ('Slovakia', 'Austria'), ('Slovakia', 'Czech Republic'),
        ('Slovakia', 'Hungary'), ('Slovakia', 'Poland'),
        ('Slovakia', 'Ukraine'), ('Slovenia', 'Austria'),
        ('Slovenia', 'Croatia'), ('Slovenia', 'Hungary'), ('Slovenia', 'Italy'),
        ('Somalia', 'Djibouti'), ('South Africa', 'Botswana'),
        ('South Africa', 'Mozambique'), ('South Africa', 'Namibia'),
        ('South Africa', 'Zimbabwe'), ('South Korea', 'North Korea'),
        ('Sudan', 'Democratic Republic of the Congo'), ('Sudan', 'Egypt'),
        ('Sudan', 'Eritrea'), ('Suriname', 'Guyana'),
        ('Swaziland', 'Mozambique'), ('Swaziland', 'South Africa'),
        ('Switzerland', 'Austria'), ('Switzerland', 'Germany'),
        ('Switzerland', 'Italy'), ('Syria', 'Iraq'), ('Syria', 'Israel'),
        ('Syria', 'Lebanon'), ('Tajikistan', 'Afghanistan'),
        ('Tajikistan', 'China'), ('Tajikistan', 'Kyrgyzstan'),
        ('Tajikistan', 'Uzbekistan'), ('Tanzania', 'Burundi'),
        ('Tanzania', 'Democratic Republic of the Congo'), ('Tanzania', 'Kenya'),
        ('Tanzania', 'Malawi'), ('Tanzania', 'Mozambique'),
        ('Tanzania', 'Rwanda'), ('Tanzania', 'Uganda'), ('Tanzania', 'Zambia'),
        ('Thailand', 'Burma'), ('The Gambia', 'Senegal'),
        ('Timor-Leste', 'Indonesia'), ('Turkey', 'Armenia'),
        ('Turkey', 'Bulgaria'), ('Turkey', 'Georgia'), ('Turkey', 'Greece'),
        ('Turkey', 'Iran'), ('Turkey', 'Iraq'), ('Turkey', 'Syria'),
        ('Turkmenistan', 'Afghanistan'), ('Turkmenistan', 'Iran'),
        ('Turkmenistan', 'Kazakhstan'), ('Turkmenistan', 'Uzbekistan'),
        ('Uganda', 'Democratic Republic of the Congo'), ('Uganda', 'Sudan'),
        ('Ukraine', 'Belarus'), ('Ukraine', 'Hungary'), ('Ukraine', 'Moldova'),
        ('Ukraine', 'Poland'), ('Ukraine', 'Romania'), ('Ukraine', 'Russia'),
        ('United States', 'Mexico'), ('Uruguay', 'Brazil'),
        ('Uzbekistan', 'Kazakhstan'), ('Uzbekistan', 'Kyrgyzstan'),
        ('Vatican City', 'Italy'), ('Venezuela', 'Guyana'),
        ('West Bank', 'Israel'), ('Western Sahara', 'Algeria'),
        ('Western Sahara', 'Morocco'), ('Zambia', 'Malawi'),
        ('Zambia', 'Zimbabwe'), ('Zimbabwe', 'Botswana')
        ]
    gps_coordinates = {
        "Cote d'Ivoire": [[8, 'N'], [5, 'W']],
        'Afghanistan': [[33, 'N'], [65, 'E']],
        'Akrotiri': [[34, 'N'], [32, 'E']],
        'Albania': [[41, 'N'], [20, 'E']],
        'Algeria': [[28, 'N'], [3, 'E']],
        'American Samoa': [[14, 'S'], [170, 'W']],
        'Andorra': [[42, 'N'], [1, 'E']],
        'Angola': [[12, 'S'], [18, 'E']],
        'Anguilla': [[18, 'N'], [63, 'W']],
        'Antarctica': [[90, 'S'], [0, 'E']],
        'Antigua and Barbuda': [[17, 'N'], [61, 'W']],
        'Argentina': [[34, 'S'], [64, 'W']],
        'Armenia': [[40, 'N'], [45, 'E']],
        'Aruba': [[12, 'N'], [69, 'W']],
        'Ashmore and Cartier Islands': [[12, 'S'], [123, 'E']],
        'Australia': [[27, 'S'], [133, 'E']],
        'Austria': [[47, 'N'], [13, 'E']],
        'Azerbaijan': [[40, 'N'], [47, 'E']],
        'Bahamas, The': [[24, 'N'], [76, 'W']],
        'Bahrain': [[26, 'N'], [50, 'E']],
        'Bangladesh': [[24, 'N'], [90, 'E']],
        'Barbados': [[13, 'N'], [59, 'W']],
        'Belarus': [[53, 'N'], [28, 'E']],
        'Belgium': [[50, 'N'], [4, 'E']],
        'Belize': [[17, 'N'], [88, 'W']],
        'Benin': [[9, 'N'], [2, 'E']],
        'Bermuda': [[32, 'N'], [64, 'W']],
        'Bhutan': [[27, 'N'], [90, 'E']],
        'Bolivia': [[17, 'S'], [65, 'W']],
        'Bosnia and Herzegovina': [[44, 'N'], [18, 'E']],
        'Botswana': [[22, 'S'], [24, 'E']],
        'Bouvet Island': [[54, 'S'], [3, 'E']],
        'Brazil': [[10, 'S'], [55, 'W']],
        'British Indian Ocean Territory': [[6, 'S'], [71, 'E']],
        'British Virgin Islands': [[18, 'N'], [64, 'W']],
        'Brunei': [[4, 'N'], [114, 'E']],
        'Bulgaria': [[43, 'N'], [25, 'E']],
        'Burkina Faso': [[13, 'N'], [2, 'W']],
        'Burma': [[22, 'N'], [98, 'E']],
        'Burundi': [[3, 'S'], [30, 'E']],
        'Cambodia': [[13, 'N'], [105, 'E']],
        'Cameroon': [[6, 'N'], [12, 'E']],
        'Canada': [[60, 'N'], [95, 'W']],
        'Cape Verde': [[16, 'N'], [24, 'W']],
        'Cayman Islands': [[19, 'N'], [80, 'W']],
        'Central African Republic': [[7, 'N'], [21, 'E']],
        'Chad': [[15, 'N'], [19, 'E']],
        'Chile': [[30, 'S'], [71, 'W']],
        'China': [[35, 'N'], [105, 'E']],
        'Christmas Island': [[10, 'S'], [105, 'E']],
        'Clipperton Island': [[10, 'N'], [109, 'W']],
        'Cocos (Keeling) Islands': [[12, 'S'], [96, 'E']],
        'Colombia': [[4, 'N'], [72, 'W']],
        'Comoros': [[12, 'S'], [44, 'E']],
        'Cook Islands': [[21, 'S'], [159, 'W']],
        'Coral Sea Islands': [[18, 'S'], [152, 'E']],
        'Costa Rica': [[10, 'N'], [84, 'W']],
        'Croatia': [[45, 'N'], [15, 'E']],
        'Cuba': [[21, 'N'], [80, 'W']],
        'Cyprus': [[35, 'N'], [33, 'E']],
        'Czech Republic': [[49, 'N'], [15, 'E']],
        'Democratic Republic of the Congo': [[0, 'N'], [25, 'E']],
        'Denmark': [[56, 'N'], [10, 'E']],
        'Dhekelia': [[34, 'N'], [33, 'E']],
        'Djibouti': [[11, 'N'], [43, 'E']],
        'Dominica': [[15, 'N'], [61, 'W']],
        'Dominican Republic': [[19, 'N'], [70, 'W']],
        'Ecuador': [[2, 'S'], [77, 'W']],
        'Egypt': [[27, 'N'], [30, 'E']],
        'El Salvador': [[13, 'N'], [88, 'W']],
        'Equatorial Guinea': [[2, 'N'], [10, 'E']],
        'Eritrea': [[15, 'N'], [39, 'E']],
        'Estonia': [[59, 'N'], [26, 'E']],
        'Ethiopia': [[8, 'N'], [38, 'E']],
        'Falkland Islands (Islas Malvinas)': [[51, 'S'], [59, 'W']],
        'Faroe Islands': [[62, 'N'], [7, 'W']],
        'Fiji': [[18, 'S'], [175, 'E']],
        'Finland': [[64, 'N'], [26, 'E']],
        'France': [[46, 'N'], [2, 'E']],
        'French Polynesia': [[15, 'S'], [140, 'W']],
        'Gabon': [[1, 'S'], [11, 'E']],
        'Gaza Strip': [[31, 'N'], [34, 'E']],
        'Georgia': [[42, 'N'], [43, 'E']],
        'Germany': [[51, 'N'], [9, 'E']],
        'Ghana': [[8, 'N'], [2, 'W']],
        'Gibraltar': [[36, 'N'], [5, 'W']],
        'Greece': [[39, 'N'], [22, 'E']],
        'Greenland': [[72, 'N'], [40, 'W']],
        'Grenada': [[12, 'N'], [61, 'W']],
        'Guam': [[13, 'N'], [144, 'E']],
        'Guatemala': [[15, 'N'], [90, 'W']],
        'Guernsey': [[49, 'N'], [2, 'W']],
        'Guinea': [[11, 'N'], [10, 'W']],
        'Guinea-Bissau': [[12, 'N'], [15, 'W']],
        'Guyana': [[5, 'N'], [59, 'W']],
        'Haiti': [[19, 'N'], [72, 'W']],
        'Heard Island and McDonald Islands': [[53, 'S'], [72, 'E']],
        'Honduras': [[15, 'N'], [86, 'W']],
        'Hong Kong': [[22, 'N'], [114, 'E']],
        'Hungary': [[47, 'N'], [20, 'E']],
        'Iceland': [[65, 'N'], [18, 'W']],
        'India': [[20, 'N'], [77, 'E']],
        'Indonesia': [[5, 'S'], [120, 'E']],
        'Iran': [[32, 'N'], [53, 'E']],
        'Iraq': [[33, 'N'], [44, 'E']],
        'Ireland': [[53, 'N'], [8, 'W']],
        'Isle of Man': [[54, 'N'], [4, 'W']],
        'Israel': [[31, 'N'], [34, 'E']],
        'Italy': [[42, 'N'], [12, 'E']],
        'Jamaica': [[18, 'N'], [77, 'W']],
        'Jan Mayen': [[71, 'N'], [8, 'W']],
        'Japan': [[36, 'N'], [138, 'E']],
        'Jersey': [[49, 'N'], [2, 'W']],
        'Jordan': [[31, 'N'], [36, 'E']],
        'Kazakhstan': [[48, 'N'], [68, 'E']],
        'Kenya': [[1, 'N'], [38, 'E']],
        'Kiribati': [[1, 'N'], [173, 'E']],
        'Kosovo': [[42, 'N'], [21, 'E']],
        'Kuwait': [[29, 'N'], [45, 'E']],
        'Kyrgyzstan': [[41, 'N'], [75, 'E']],
        'Laos': [[18, 'N'], [105, 'E']],
        'Latvia': [[57, 'N'], [25, 'E']],
        'Lebanon': [[33, 'N'], [35, 'E']],
        'Lesotho': [[29, 'S'], [28, 'E']],
        'Liberia': [[6, 'N'], [9, 'W']],
        'Libya': [[25, 'N'], [17, 'E']],
        'Liechtenstein': [[47, 'N'], [9, 'E']],
        'Lithuania': [[56, 'N'], [24, 'E']],
        'Luxembourg': [[49, 'N'], [6, 'E']],
        'Macau': [[22, 'N'], [113, 'E']],
        'Macedonia': [[41, 'N'], [22, 'E']],
        'Madagascar': [[20, 'S'], [47, 'E']],
        'Malawi': [[13, 'S'], [34, 'E']],
        'Malaysia': [[2, 'N'], [112, 'E']],
        'Maldives': [[3, 'N'], [73, 'E']],
        'Mali': [[17, 'N'], [4, 'W']],
        'Malta': [[35, 'N'], [14, 'E']],
        'Marshall Islands': [[9, 'N'], [168, 'E']],
        'Mauritania': [[20, 'N'], [12, 'W']],
        'Mauritius': [[20, 'S'], [57, 'E']],
        'Mayotte': [[12, 'S'], [45, 'E']],
        'Mexico': [[23, 'N'], [102, 'W']],
        'Micronesia, Federated States of': [[6, 'N'], [158, 'E']],
        'Moldova': [[47, 'N'], [29, 'E']],
        'Monaco': [[43, 'N'], [7, 'E']],
        'Mongolia': [[46, 'N'], [105, 'E']],
        'Montenegro': [[42, 'N'], [19, 'E']],
        'Montserrat': [[16, 'N'], [62, 'W']],
        'Morocco': [[32, 'N'], [5, 'W']],
        'Mozambique': [[18, 'S'], [35, 'E']],
        'Namibia': [[22, 'S'], [17, 'E']],
        'Nauru': [[0, 'S'], [166, 'E']],
        'Navassa Island': [[18, 'N'], [75, 'W']],
        'Nepal': [[28, 'N'], [84, 'E']],
        'Netherlands Antilles': [[12, 'N'], [69, 'W']],
        'Netherlands': [[52, 'N'], [5, 'E']],
        'New Caledonia': [[21, 'S'], [165, 'E']],
        'New Zealand': [[41, 'S'], [174, 'E']],
        'Nicaragua': [[13, 'N'], [85, 'W']],
        'Niger': [[16, 'N'], [8, 'E']],
        'Nigeria': [[10, 'N'], [8, 'E']],
        'Niue': [[19, 'S'], [169, 'W']],
        'Norfolk Island': [[29, 'S'], [167, 'E']],
        'North Korea': [[40, 'N'], [127, 'E']],
        'Northern Mariana Islands': [[15, 'N'], [145, 'E']],
        'Norway': [[62, 'N'], [10, 'E']],
        'Oman': [[21, 'N'], [57, 'E']],
        'Pakistan': [[30, 'N'], [70, 'E']],
        'Palau': [[7, 'N'], [134, 'E']],
        'Panama': [[9, 'N'], [80, 'W']],
        'Papua New Guinea': [[6, 'S'], [147, 'E']],
        'Paracel Islands': [[16, 'N'], [112, 'E']],
        'Paraguay': [[23, 'S'], [58, 'W']],
        'Peru': [[10, 'S'], [76, 'W']],
        'Philippines': [[13, 'N'], [122, 'E']],
        'Pitcairn Islands': [[25, 'S'], [130, 'W']],
        'Poland': [[52, 'N'], [20, 'E']],
        'Portugal': [[39, 'N'], [8, 'W']],
        'Puerto Rico': [[18, 'N'], [66, 'W']],
        'Qatar': [[25, 'N'], [51, 'E']],
        'Republic of the Congo': [[1, 'S'], [15, 'E']],
        'Romania': [[46, 'N'], [25, 'E']],
        'Russia': [[60, 'N'], [100, 'E']],
        'Rwanda': [[2, 'S'], [30, 'E']],
        'Saint Barthelemy': [[17, 'N'], [62, 'W']],
        'Saint Helena': [[15, 'S'], [5, 'W']],
        'Saint Kitts and Nevis': [[17, 'N'], [62, 'W']],
        'Saint Lucia': [[13, 'N'], [60, 'W']],
        'Saint Martin': [[18, 'N'], [63, 'W']],
        'Saint Pierre and Miquelon': [[46, 'N'], [56, 'W']],
        'Saint Vincent and the Grenadines': [[13, 'N'], [61, 'W']],
        'Samoa': [[13, 'S'], [172, 'W']],
        'San Marino': [[43, 'N'], [12, 'E']],
        'Sao Tome and Principe': [[1, 'N'], [7, 'E']],
        'Saudi Arabia': [[25, 'N'], [45, 'E']],
        'Senegal': [[14, 'N'], [14, 'W']],
        'Serbia': [[44, 'N'], [21, 'E']],
        'Seychelles': [[4, 'S'], [55, 'E']],
        'Sierra Leone': [[8, 'N'], [11, 'W']],
        'Singapore': [[1, 'N'], [103, 'E']],
        'Slovakia': [[48, 'N'], [19, 'E']],
        'Slovenia': [[46, 'N'], [14, 'E']],
        'Solomon Islands': [[8, 'S'], [159, 'E']],
        'Somalia': [[10, 'N'], [49, 'E']],
        'South Africa': [[29, 'S'], [24, 'E']],
        'South Georgia and the South Sandwich Islands': [[54, 'S'], [37, 'W']],
        'South Korea': [[37, 'N'], [127, 'E']],
        'Spain': [[40, 'N'], [4, 'W']],
        'Spratly Islands': [[8, 'N'], [111, 'E']],
        'Sri Lanka': [[7, 'N'], [81, 'E']],
        'Sudan': [[15, 'N'], [30, 'E']],
        'Suriname': [[4, 'N'], [56, 'W']],
        'Svalbard': [[78, 'N'], [20, 'E']],
        'Swaziland': [[26, 'S'], [31, 'E']],
        'Sweden': [[62, 'N'], [15, 'E']],
        'Switzerland': [[47, 'N'], [8, 'E']],
        'Syria': [[35, 'N'], [38, 'E']],
        'Taiwan': [[23, 'N'], [121, 'E']],
        'Tajikistan': [[39, 'N'], [71, 'E']],
        'Tanzania': [[6, 'S'], [35, 'E']],
        'Thailand': [[15, 'N'], [100, 'E']],
        'The Gambia': [[13, 'N'], [16, 'W']],
        'Timor-Leste': [[8, 'S'], [125, 'E']],
        'Togo': [[8, 'N'], [1, 'E']],
        'Tokelau': [[9, 'S'], [172, 'W']],
        'Tonga': [[20, 'S'], [175, 'W']],
        'Trinidad and Tobago': [[11, 'N'], [61, 'W']],
        'Tunisia': [[34, 'N'], [9, 'E']],
        'Turkey': [[39, 'N'], [35, 'E']],
        'Turkmenistan': [[40, 'N'], [60, 'E']],
        'Turks and Caicos Islands': [[21, 'N'], [71, 'W']],
        'Tuvalu': [[8, 'S'], [178, 'E']],
        'Uganda': [[1, 'N'], [32, 'E']],
        'Ukraine': [[49, 'N'], [32, 'E']],
        'United Arab Emirates': [[24, 'N'], [54, 'E']],
        'United Kingdom': [[54, 'N'], [2, 'W']],
        'United States': [[38, 'N'], [97, 'W']],
        'Uruguay': [[33, 'S'], [56, 'W']],
        'Uzbekistan': [[41, 'N'], [64, 'E']],
        'Vanuatu': [[16, 'S'], [167, 'E']],
        'Vatican City': [[41, 'N'], [12, 'E']],
        'Venezuela': [[8, 'N'], [66, 'W']],
        'Vietnam': [[16, 'N'], [107, 'E']],
        'Virgin Islands': [[18, 'N'], [64, 'W']],
        'Wake Island': [[19, 'N'], [166, 'E']],
        'Wallis and Futuna': [[13, 'S'], [176, 'W']],
        'West Bank': [[32, 'N'], [35, 'E']],
        'Western Sahara': [[24, 'N'], [13, 'W']],
        'Yemen': [[15, 'N'], [48, 'E']],
        'Zambia': [[15, 'S'], [30, 'E']],
        'Zimbabwe': [[20, 'S'], [30, 'E']]
        }
    g = Graph()
    g.add_edges(edges)
    g.add_vertices(gps_coordinates)
    g.gps_coordinates = gps_coordinates
    g.name("World Map")
    return g
