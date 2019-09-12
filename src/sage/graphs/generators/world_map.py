# -*- coding: utf-8 -*-
r"""
Map graphs

The methods defined here appear in :mod:`sage.graphs.graph_generators`.
"""

###########################################################################
#
#           Copyright (C) 2006 Robert L. Miller <rlmillster@gmail.com>
#                              and Emily A. Kirkman
#           Copyright (C) 2009 Michael C. Yurko <myurko@gmail.com>
#
# Distributed  under  the  terms  of  the  GNU  General  Public  License (GPL)
#                         http://www.gnu.org/licenses/
###########################################################################

# import from Sage library
from sage.graphs.graph import Graph

def AfricaMap(continental=False, year=2018):
    """
    Return African states as a graph of common border.

    "African state" here is defined as an independent
    state having the capital city in Africa. The graph
    has an edge between those countries that have common
    *land* border.

    INPUT:

    - ``continental``, a Boolean -- if set, only return states in
      the continental Africa
    - ``year`` -- reserved for future use

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

        sage: Africa.plot()
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
     'Cameroon': ['Central Africa', 'Chad', 'Equatorial Guinea', 'Gabon', 'Nigeria'],
     'Central Africa': ['Chad', 'South Sudan', 'Sudan'],
     'Chad': ['Libya', 'Niger', 'Nigeria', 'Sudan'],
     'Republic of the Congo': ['Gabon', 'Cameroon', 'Central Africa', 'Angola',
                               'Democratic Republic of the Congo'],
     'Democratic Republic of the Congo': ['Zambia', 'South Sudan', 'Tanzania', 'Burundi',
                                          'Rwanda', 'Uganda', 'Central Africa', 'Angola'],
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

    no_land_border = ['Cape Verde', 'Seychelles', 'Mauritius', u'São Tomé and Príncipe', 'Madagascar', 'Comoros']

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

    "European state" here is defined as an independent
    state having the capital city in Europe. The graph
    has an edge between those countries that have common
    *land* border.

    INPUT:

    - ``continental``, a Boolean -- if set, only return states in
      the continental Europe
    - ``year`` -- reserved for future use

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
     'Poland': ['Slovakia', 'Czech Republic', 'Lithuania', 'Russia', 'Ukraine', 'Germany'],
     'Germany': ['Czech Republic', 'Netherlands', 'Switzerland', 'Luxembourg', 'Denmark'],
     'Croatia': ['Bosnia and Herzegovina', 'Serbia', 'Hungary', 'Montenegro', 'Slovenia'],
     'Austria': ['Czech Republic', 'Germany', 'Switzerland', 'Slovenia', 'Liechtenstein'],
     'France': ['Germany', 'Italy', 'Switzerland', 'Monaco', 'Luxembourg', 'Andorra'],
     'Hungary': ['Slovakia', 'Serbia', 'Romania', 'Ukraine', 'Slovenia', 'Austria'],
     'Italy': ['Switzerland', 'Vatican City', 'San Marino', 'Slovenia', 'Austria'],
     'Belarus': ['Poland', 'Latvia', 'Lithuania', 'Russia', 'Ukraine'],
     'Montenegro': ['Bosnia and Herzegovina', 'Serbia', 'Albania'],
     'Belgium': ['Germany', 'Netherlands', 'Luxembourg', 'France'],
     'Russia': ['Finland', 'Lithuania', 'Estonia', 'Ukraine'],
     'Romania': ['Serbia', 'Moldova', 'Bulgaria', 'Ukraine'],
     'Latvia': ['Lithuania', 'Russia', 'Estonia'],
     'Slovakia': ['Czech Republic', 'Ukraine', 'Austria'], 'Switzerland': ['Liechtenstein'],
     'Spain': ['Portugal', 'Andorra', 'France'], 'Norway': ['Finland', 'Sweden', 'Russia'],
     'Ireland': ['United Kingdom'], 'Serbia': ['Bosnia and Herzegovina', 'Bulgaria'],
     'Greece': ['Macedonia', 'Bulgaria', 'Albania'], 'Ukraine': ['Moldova'],
     'Macedonia': ['Serbia', 'Bulgaria', 'Albania'], 'Sweden': ['Finland']
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

    The graph has an edge between those states that have
    common *land* border line or point. Hence for example
    Colorado and Arizona are marked as neighbors, but
    Michigan and Minnesota are not.

    INPUT:

    - ``continental``, a Boolean -- if set, exclude Alaska
      and Hawaii

    EXAMPLES:

    How many states are neighbor's neighbor for Pennsylvania::

        sage: USA = graphs.USAMap()
        sage: len([n2 for n2 in USA if USA.distance('Pennsylvania', n2) == 2])
        7

    Diameter for continental USA::

        sage: USAcont = graphs.USAMap(continental=True)
        sage: USAcont.diameter()
        11
    """
    states = {
    "Alabama": ["Florida", "Georgia", "Mississippi", "Tennessee"],
    "Arizona": ["California", "Colorado", "Nevada", "New Mexico", "Utah"],
    "Arkansas": ["Louisiana", "Mississippi", "Missouri", "Oklahoma", "Tennessee", "Texas"],
    "California": ["Arizona", "Nevada", "Oregon"],
    "Colorado": ["Arizona", "Kansas", "Nebraska", "New Mexico", "Oklahoma", "Utah", "Wyoming"],
    "Connecticut": ["Massachusetts", "New York", "Rhode Island"],
    "Delaware": ["Maryland", "New Jersey", "Pennsylvania"],
    "Florida": ["Alabama", "Georgia"],
    "Georgia": ["Alabama", "Florida", "North Carolina", "South Carolina", "Tennessee"],
    "Idaho": ["Montana", "Nevada", "Oregon", "Utah", "Washington", "Wyoming"],
    "Illinois": ["Indiana", "Iowa", "Michigan", "Kentucky", "Missouri", "Wisconsin"],
    "Indiana": ["Illinois", "Kentucky", "Michigan", "Ohio"],
    "Iowa": ["Illinois", "Minnesota", "Missouri", "Nebraska", "South Dakota", "Wisconsin"],
    "Kansas": ["Colorado", "Missouri", "Nebraska", "Oklahoma"],
    "Kentucky": ["Illinois", "Indiana", "Missouri", "Ohio", "Tennessee", "Virginia", "West Virginia"],
    "Louisiana": ["Arkansas", "Mississippi", "Texas"],
    "Maine": ["New Hampshire"],
    "Maryland": ["Delaware", "Pennsylvania", "Virginia", "West Virginia"],
    "Massachusetts": ["Connecticut", "New Hampshire", "New York", "Rhode Island", "Vermont"],
    "Michigan": ["Illinois", "Indiana", "Ohio", "Wisconsin"],
    "Minnesota": ["Iowa", "North Dakota", "South Dakota", "Wisconsin"],
    "Mississippi": ["Alabama", "Arkansas", "Louisiana", "Tennessee"],
    "Missouri": ["Arkansas", "Illinois", "Iowa", "Kansas", "Kentucky", "Nebraska", "Oklahoma", "Tennessee"],
    "Montana": ["Idaho", "North Dakota", "South Dakota", "Wyoming"],
    "Nebraska": ["Colorado", "Iowa", "Kansas", "Missouri", "South Dakota", "Wyoming"],
    "Nevada": ["Arizona", "California", "Idaho", "Oregon", "Utah"],
    "New Hampshire": ["Maine", "Massachusetts", "Vermont"],
    "New Jersey": ["Delaware", "New York", "Pennsylvania"],
    "New Mexico": ["Arizona", "Colorado", "Oklahoma", "Texas", "Utah"],
    "New York": ["Connecticut", "Massachusetts", "New Jersey", "Pennsylvania", "Vermont"],
    "North Carolina": ["Georgia", "South Carolina", "Tennessee", "Virginia"],
    "North Dakota": ["Minnesota", "Montana", "South Dakota"],
    "Ohio": ["Indiana", "Kentucky", "Michigan", "Pennsylvania", "West Virginia"],
    "Oklahoma": ["Arkansas", "Colorado", "Kansas", "Missouri", "New Mexico", "Texas"],
    "Oregon": ["California", "Idaho", "Nevada", "Washington"],
    "Pennsylvania": ["Delaware", "Maryland", "New Jersey", "New York", "Ohio", "West Virginia"],
    "Rhode Island": ["Connecticut", "Massachusetts"],
    "South Carolina": ["Georgia", "North Carolina"],
    "South Dakota": ["Iowa", "Minnesota", "Montana", "Nebraska", "North Dakota", "Wyoming"],
    "Tennessee": ["Alabama", "Arkansas", "Georgia", "Kentucky", "Mississippi", "Missouri", "North Carolina", "Virginia"],
    "Texas": ["Arkansas", "Louisiana", "New Mexico", "Oklahoma"],
    "Utah": ["Arizona", "Colorado", "Idaho", "Nevada", "New Mexico", "Wyoming"],
    "Vermont": ["Massachusetts", "New Hampshire", "New York"],
    "Virginia": ["Kentucky", "Maryland", "North Carolina", "Tennessee", "West Virginia"],
    "Washington": ["Idaho", "Oregon"],
    "West Virginia": ["Kentucky", "Maryland", "Ohio", "Pennsylvania", "Virginia"],
    "Wisconsin": ["Illinois", "Iowa", "Michigan", "Minnesota"],
    "Wyoming": ["Colorado", "Idaho", "Montana", "Nebraska", "South Dakota", "Utah"]
    }
    if not continental:
        states['Alaska'] = []
        states['Hawaii'] = []
        G = Graph(states, format='dict_of_lists')
        G.name(new="USA Map")
        return G

    G = Graph(states, format='dict_of_lists')
    G.name(new="Continental USA Map")
    return G

def WorldMap():
    """
    Returns the Graph of all the countries, in which two countries are adjacent
    in the graph if they have a common boundary.

    This graph has been built from the data available
    in The CIA World Factbook [CIA]_ (2009-08-21).

    The returned graph ``G`` has a member ``G.gps_coordinates``
    equal to a dictionary containing the GPS coordinates
    of each country's capital city.

    EXAMPLES::

        sage: g = graphs.WorldMap()
        sage: g.has_edge("France", "Italy")
        True
        sage: g.gps_coordinates["Bolivia"]
        [[17, 'S'], [65, 'W']]
        sage: sorted(g.connected_component_containing_vertex('Ireland'))
        ['Ireland', 'United Kingdom']

    TESTS::

        sage: 'Iceland' in graphs.WorldMap()  # Trac 24488
        True

    REFERENCE:

    [CIA]_
    """
    edges = [
        ('Afghanistan', 'China', None), ('Afghanistan', 'Iran', None),
        ('Afghanistan', 'Uzbekistan', None), ('Albania', 'Greece', None),
        ('Albania', 'Kosovo', None), ('Albania', 'Macedonia', None),
        ('Albania', 'Montenegro', None), ('Algeria', 'Morocco', None),
        ('Algeria', 'Tunisia', None), ('Andorra', 'Spain', None),
        ('Angola', 'Democratic Republic of the Congo', None), ('Angola', 'Namibia', None),
        ('Angola', 'Zambia', None), ('Argentina', 'Bolivia', None),
        ('Argentina', 'Brazil', None), ('Argentina', 'Chile', None),
        ('Argentina', 'Paraguay', None), ('Argentina', 'Uruguay', None),
        ('Armenia', 'Georgia', None), ('Armenia', 'Iran', None),
        ('Austria', 'Germany', None), ('Azerbaijan', 'Armenia', None),
        ('Azerbaijan', 'Georgia', None), ('Azerbaijan', 'Iran', None),
        ('Azerbaijan', 'Russia', None), ('Azerbaijan', 'Turkey', None),
        ('Bangladesh', 'Burma', None), ('Belgium', 'Germany', None),
        ('Belgium', 'Netherlands', None), ('Belize', 'Mexico', None),
        ('Benin', 'Burkina Faso', None), ('Benin', 'Niger', None),
        ('Benin', 'Nigeria', None), ('Benin', 'Togo', None),
        ('Bolivia', 'Brazil', None), ('Bolivia', 'Chile', None),
        ('Bolivia', 'Paraguay', None), ('Bolivia', 'Peru', None),
        ('Bosnia and Herzegovina', 'Croatia', None), ('Bosnia and Herzegovina', 'Montenegro', None),
        ('Bosnia and Herzegovina', 'Serbia', None), ('Brazil', 'Colombia', None),
        ('Brazil', 'Guyana', None), ('Brazil', 'Suriname', None),
        ('Brazil', 'Venezuela', None), ('Bulgaria', 'Greece', None),
        ('Bulgaria', 'Macedonia', None), ('Bulgaria', 'Romania', None),
        ('Bulgaria', 'Serbia', None), ('Burkina Faso', 'Mali', None),
        ('Burkina Faso', 'Niger', None), ('Burkina Faso', 'Togo', None),
        ('Burundi', 'Democratic Republic of the Congo', None), ('Cambodia', 'Laos', None),
        ('Cambodia', 'Thailand', None), ('Cambodia', 'Vietnam', None),
        ('Cameroon', 'Central African Republic', None), ('Cameroon', 'Chad', None),
        ('Cameroon', 'Equatorial Guinea', None), ('Cameroon', 'Nigeria', None),
        ('Cameroon', 'Republic of the Congo', None), ('Canada', 'United States', None),
        ('Central African Republic', 'Chad', None), ('Central African Republic', 'Democratic Republic of the Congo', None),
        ('Central African Republic', 'Sudan', None), ('Chad', 'Niger', None),
        ('Chad', 'Nigeria', None), ('Chad', 'Sudan', None),
        ('China', 'Bhutan', None), ('China', 'Burma', None),
        ('China', 'Hong Kong', None), ('China', 'Kazakhstan', None),
        ('China', 'Kyrgyzstan', None), ('China', 'Mongolia', None),
        ('China', 'Nepal', None), ('China', 'North Korea', None),
        ('China', 'Russia', None), ('China', 'Vietnam', None),
        ('Colombia', 'Venezuela', None), ('Costa Rica', 'Nicaragua', None),
        ("Cote d'Ivoire", 'Burkina Faso', None), ("Cote d'Ivoire", 'Guinea', None),
        ("Cote d'Ivoire", 'Mali', None), ('Cyprus', 'Akrotiri', None),
        ('Cyprus', 'Dhekelia', None), ('Czech Republic', 'Austria', None),
        ('Czech Republic', 'Germany', None), ('Czech Republic', 'Poland', None),
        ('Democratic Republic of the Congo', 'Zambia', None), ('Denmark', 'Germany', None),
        ('Djibouti', 'Eritrea', None), ('Dominican Republic', 'Haiti', None),
        ('Ecuador', 'Colombia', None), ('El Salvador', 'Honduras', None),
        ('Ethiopia', 'Djibouti', None), ('Ethiopia', 'Eritrea', None),
        ('Ethiopia', 'Kenya', None), ('Ethiopia', 'Somalia', None),
        ('Ethiopia', 'Sudan', None), ('Finland', 'Russia', None),
        ('Finland', 'Sweden', None), ('France', 'Andorra', None),
        ('France', 'Belgium', None), ('France', 'Brazil', None),
        ('France', 'Germany', None), ('France', 'Italy', None),
        ('France', 'Luxembourg', None), ('France', 'Spain', None),
        ('France', 'Suriname', None), ('France', 'Switzerland', None),
        ('Gabon', 'Cameroon', None), ('Gabon', 'Equatorial Guinea', None),
        ('Gabon', 'Republic of the Congo', None), ('Gaza Strip', 'Egypt', None),
        ('Gaza Strip', 'Israel', None), ('Ghana', 'Burkina Faso', None),
        ('Ghana', "Cote d'Ivoire", None), ('Ghana', 'Togo', None),
        ('Gibraltar', 'Spain', None), ('Guatemala', 'Belize', None),
        ('Guatemala', 'El Salvador', None), ('Guatemala', 'Honduras', None),
        ('Guatemala', 'Mexico', None), ('Guinea', 'Sierra Leone', None),
        ('Guinea-Bissau', 'Guinea', None), ('Guinea-Bissau', 'Senegal', None),
        ('Honduras', 'Nicaragua', None), ('Hungary', 'Austria', None),
        ('Hungary', 'Croatia', None), ('Hungary', 'Serbia', None),
        ('India', 'Bangladesh', None), ('India', 'Bhutan', None),
        ('India', 'Burma', None), ('India', 'China', None),
        ('India', 'Nepal', None), ('Indonesia', 'Papua New Guinea', None),
        ('Iran', 'Iraq', None), ('Ireland', 'United Kingdom', None),
        ('Israel', 'Egypt', None), ('Italy', 'Austria', None),
        ('Jordan', 'Iraq', None), ('Jordan', 'Israel', None),
        ('Jordan', 'Syria', None), ('Jordan', 'West Bank', None),
        ('Kazakhstan', 'Kyrgyzstan', None), ('Kenya', 'Somalia', None),
        ('Kenya', 'Sudan', None), ('Kenya', 'Uganda', None),
        ('Kosovo', 'Macedonia', None), ('Kosovo', 'Serbia', None),
        ('Kuwait', 'Iraq', None), ('Laos', 'Burma', None),
        ('Laos', 'China', None), ('Laos', 'Thailand', None),
        ('Laos', 'Vietnam', None), ('Latvia', 'Belarus', None),
        ('Latvia', 'Estonia', None), ('Lebanon', 'Israel', None),
        ('Lesotho', 'South Africa', None), ('Liberia', "Cote d'Ivoire", None),
        ('Liberia', 'Guinea', None), ('Liberia', 'Sierra Leone', None),
        ('Libya', 'Algeria', None), ('Libya', 'Chad', None),
        ('Libya', 'Egypt', None), ('Libya', 'Niger', None),
        ('Libya', 'Sudan', None), ('Libya', 'Tunisia', None),
        ('Liechtenstein', 'Austria', None), ('Liechtenstein', 'Switzerland', None),
        ('Lithuania', 'Belarus', None), ('Lithuania', 'Latvia', None),
        ('Lithuania', 'Poland', None), ('Lithuania', 'Russia', None),
        ('Luxembourg', 'Belgium', None), ('Luxembourg', 'Germany', None),
        ('Macau', 'China', None), ('Macedonia', 'Greece', None),
        ('Macedonia', 'Serbia', None), ('Malaysia', 'Brunei', None),
        ('Malaysia', 'Indonesia', None), ('Malaysia', 'Thailand', None),
        ('Mali', 'Algeria', None), ('Mali', 'Guinea', None),
        ('Mali', 'Niger', None), ('Mali', 'Senegal', None),
        ('Mauritania', 'Algeria', None), ('Mauritania', 'Mali', None),
        ('Mauritania', 'Senegal', None), ('Mauritania', 'Western Sahara', None),
        ('Monaco', 'France', None), ('Montenegro', 'Croatia', None),
        ('Montenegro', 'Kosovo', None), ('Montenegro', 'Serbia', None),
        ('Morocco', 'Spain', None), ('Mozambique', 'Malawi', None),
        ('Mozambique', 'Zambia', None), ('Mozambique', 'Zimbabwe', None),
        ('Namibia', 'Botswana', None), ('Namibia', 'Zambia', None),
        ('Netherlands', 'Germany', None), ('Niger', 'Algeria', None),
        ('Niger', 'Nigeria', None), ('Norway', 'Finland', None),
        ('Norway', 'Russia', None), ('Norway', 'Sweden', None),
        ('Oman', 'United Arab Emirates', None), ('Oman', 'Yemen', None),
        ('Pakistan', 'Afghanistan', None), ('Pakistan', 'China', None),
        ('Pakistan', 'India', None), ('Pakistan', 'Iran', None),
        ('Panama', 'Colombia', None), ('Panama', 'Costa Rica', None),
        ('Paraguay', 'Brazil', None), ('Peru', 'Brazil', None),
        ('Peru', 'Chile', None), ('Peru', 'Colombia', None),
        ('Peru', 'Ecuador', None), ('Poland', 'Belarus', None),
        ('Poland', 'Germany', None), ('Portugal', 'Spain', None),
        ('Republic of the Congo', 'Angola', None), ('Republic of the Congo', 'Central African Republic', None),
        ('Republic of the Congo', 'Democratic Republic of the Congo', None), ('Romania', 'Hungary', None),
        ('Romania', 'Moldova', None), ('Romania', 'Serbia', None),
        ('Russia', 'Belarus', None), ('Russia', 'Estonia', None),
        ('Russia', 'Georgia', None), ('Russia', 'Kazakhstan', None),
        ('Russia', 'Latvia', None), ('Russia', 'Mongolia', None),
        ('Russia', 'North Korea', None), ('Russia', 'Poland', None),
        ('Rwanda', 'Burundi', None), ('Rwanda', 'Democratic Republic of the Congo', None),
        ('Rwanda', 'Uganda', None), ('Saint Martin', 'Netherlands Antilles', None),
        ('San Marino', 'Italy', None), ('Saudi Arabia', 'Iraq', None),
        ('Saudi Arabia', 'Jordan', None), ('Saudi Arabia', 'Kuwait', None),
        ('Saudi Arabia', 'Oman', None), ('Saudi Arabia', 'Qatar', None),
        ('Saudi Arabia', 'United Arab Emirates', None), ('Saudi Arabia', 'Yemen', None),
        ('Senegal', 'Guinea', None), ('Serbia', 'Croatia', None),
        ('Slovakia', 'Austria', None), ('Slovakia', 'Czech Republic', None),
        ('Slovakia', 'Hungary', None), ('Slovakia', 'Poland', None),
        ('Slovakia', 'Ukraine', None), ('Slovenia', 'Austria', None),
        ('Slovenia', 'Croatia', None), ('Slovenia', 'Hungary', None),
        ('Slovenia', 'Italy', None), ('Somalia', 'Djibouti', None),
        ('South Africa', 'Botswana', None), ('South Africa', 'Mozambique', None),
        ('South Africa', 'Namibia', None), ('South Africa', 'Zimbabwe', None),
        ('South Korea', 'North Korea', None), ('Sudan', 'Democratic Republic of the Congo', None),
        ('Sudan', 'Egypt', None), ('Sudan', 'Eritrea', None),
        ('Suriname', 'Guyana', None), ('Swaziland', 'Mozambique', None),
        ('Swaziland', 'South Africa', None), ('Switzerland', 'Austria', None),
        ('Switzerland', 'Germany', None), ('Switzerland', 'Italy', None),
        ('Syria', 'Iraq', None), ('Syria', 'Israel', None),
        ('Syria', 'Lebanon', None), ('Tajikistan', 'Afghanistan', None),
        ('Tajikistan', 'China', None), ('Tajikistan', 'Kyrgyzstan', None),
        ('Tajikistan', 'Uzbekistan', None), ('Tanzania', 'Burundi', None),
        ('Tanzania', 'Democratic Republic of the Congo', None), ('Tanzania', 'Kenya', None),
        ('Tanzania', 'Malawi', None), ('Tanzania', 'Mozambique', None),
        ('Tanzania', 'Rwanda', None), ('Tanzania', 'Uganda', None),
        ('Tanzania', 'Zambia', None), ('Thailand', 'Burma', None),
        ('The Gambia', 'Senegal', None), ('Timor-Leste', 'Indonesia', None),
        ('Turkey', 'Armenia', None), ('Turkey', 'Bulgaria', None),
        ('Turkey', 'Georgia', None), ('Turkey', 'Greece', None),
        ('Turkey', 'Iran', None), ('Turkey', 'Iraq', None),
        ('Turkey', 'Syria', None), ('Turkmenistan', 'Afghanistan', None),
        ('Turkmenistan', 'Iran', None), ('Turkmenistan', 'Kazakhstan', None),
        ('Turkmenistan', 'Uzbekistan', None), ('Uganda', 'Democratic Republic of the Congo', None),
        ('Uganda', 'Sudan', None), ('Ukraine', 'Belarus', None),
        ('Ukraine', 'Hungary', None), ('Ukraine', 'Moldova', None),
        ('Ukraine', 'Poland', None), ('Ukraine', 'Romania', None),
        ('Ukraine', 'Russia', None), ('United States', 'Mexico', None),
        ('Uruguay', 'Brazil', None), ('Uzbekistan', 'Kazakhstan', None),
        ('Uzbekistan', 'Kyrgyzstan', None), ('Vatican City', 'Italy', None),
        ('Venezuela', 'Guyana', None), ('West Bank', 'Israel', None),
        ('Western Sahara', 'Algeria', None), ('Western Sahara', 'Morocco', None),
        ('Zambia', 'Malawi', None), ('Zambia', 'Zimbabwe', None),
        ('Zimbabwe', 'Botswana', None)
        ]
    gps_coordinates = {
        'Canada': [[60, 'N'], [95, 'W']],
        'Saint Martin': [[18, 'N'], [63, 'W']],
        'Sao Tome and Principe': [[1, 'N'], [7, 'E']],
        'Turkmenistan': [[40, 'N'], [60, 'E']],
        'Saint Helena': [[15, 'S'], [5, 'W']],
        'Lithuania': [[56, 'N'], [24, 'E']],
        'Cambodia': [[13, 'N'], [105, 'E']],
        'Saint Kitts and Nevis': [[17, 'N'], [62, 'W']],
        'Ethiopia': [[8, 'N'], [38, 'E']],
        'The Gambia': [[13, 'N'], [16, 'W']],
        'Aruba': [[12, 'N'], [69, 'W']],
        'Swaziland': [[26, 'S'], [31, 'E']],
        'Guinea-Bissau': [[12, 'N'], [15, 'W']],
        'Argentina': [[34, 'S'], [64, 'W']],
        'Bolivia': [[17, 'S'], [65, 'W']],
        'Bahamas, The': [[24, 'N'], [76, 'W']],
        'Spratly Islands': [[8, 'N'], [111, 'E']],
        'Ghana': [[8, 'N'], [2, 'W']],
        'Saudi Arabia': [[25, 'N'], [45, 'E']],
        'American Samoa': [[14, 'S'], [170, 'W']],
        'Cocos (Keeling) Islands': [[12, 'S'], [96, 'E']],
        'Slovenia': [[46, 'N'], [14, 'E']],
        'Guatemala': [[15, 'N'], [90, 'W']],
        'Bosnia and Herzegovina': [[44, 'N'], [18, 'E']],
        'Kuwait': [[29, 'N'], [45, 'E']],
        'Jordan': [[31, 'N'], [36, 'E']],
        'Saint Barthelemy': [[17, 'N'], [62, 'W']],
        'Ashmore and Cartier Islands': [[12, 'S'], [123, 'E']],
        'Dominica': [[15, 'N'], [61, 'W']],
        'Liberia': [[6, 'N'], [9, 'W']],
        'Maldives': [[3, 'N'], [73, 'E']],
        'Micronesia, Federated States of': [[6, 'N'], [158, 'E']],
        'Pakistan': [[30, 'N'], [70, 'E']],
        'Oman': [[21, 'N'], [57, 'E']],
        'Tanzania': [[6, 'S'], [35, 'E']],
        'Albania': [[41, 'N'], [20, 'E']],
        'Gabon': [[1, 'S'], [11, 'E']],
        'Niue': [[19, 'S'], [169, 'W']],
        'Monaco': [[43, 'N'], [7, 'E']],
        'Wallis and Futuna': [[13, 'S'], [176, 'W']],
        'New Zealand': [[41, 'S'], [174, 'E']],
        'Yemen': [[15, 'N'], [48, 'E']],
        'Jersey': [[49, 'N'], [2, 'W']],
        'Jamaica': [[18, 'N'], [77, 'W']],
        'Greenland': [[72, 'N'], [40, 'W']],
        'West Bank': [[32, 'N'], [35, 'E']],
        'Macau': [[22, 'N'], [113, 'E']],
        'Jan Mayen': [[71, 'N'], [8, 'W']],
        'United Arab Emirates': [[24, 'N'], [54, 'E']],
        'Guam': [[13, 'N'], [144, 'E']],
        'Uruguay': [[33, 'S'], [56, 'W']],
        'India': [[20, 'N'], [77, 'E']],
        'Azerbaijan': [[40, 'N'], [47, 'E']],
        'Lesotho': [[29, 'S'], [28, 'E']],
        'Saint Vincent and the Grenadines': [[13, 'N'], [61, 'W']],
        'Kenya': [[1, 'N'], [38, 'E']],
        'South Korea': [[37, 'N'], [127, 'E']],
        'Tajikistan': [[39, 'N'], [71, 'E']],
        'Turkey': [[39, 'N'], [35, 'E']],
        'Afghanistan': [[33, 'N'], [65, 'E']],
        'Paraguay': [[23, 'S'], [58, 'W']],
        'Bangladesh': [[24, 'N'], [90, 'E']],
        'Mauritania': [[20, 'N'], [12, 'W']],
        'Solomon Islands': [[8, 'S'], [159, 'E']],
        'Saint Pierre and Miquelon': [[46, 'N'], [56, 'W']],
        'Gaza Strip': [[31, 'N'], [34, 'E']],
        'San Marino': [[43, 'N'], [12, 'E']],
        'French Polynesia': [[15, 'S'], [140, 'W']],
        'France': [[46, 'N'], [2, 'E']],
        'Fiji': [[18, 'S'], [175, 'E']],
        'Rwanda': [[2, 'S'], [30, 'E']],
        'Slovakia': [[48, 'N'], [19, 'E']],
        'Somalia': [[10, 'N'], [49, 'E']],
        'Peru': [[10, 'S'], [76, 'W']],
        'Laos': [[18, 'N'], [105, 'E']],
        'Nauru': [[0, 'S'], [166, 'E']],
        'Seychelles': [[4, 'S'], [55, 'E']],
        'Norway': [[62, 'N'], [10, 'E']],
        "Cote d'Ivoire": [[8, 'N'], [5, 'W']],
        'Cook Islands': [[21, 'S'], [159, 'W']],
        'Benin': [[9, 'N'], [2, 'E']],
        'Western Sahara': [[24, 'N'], [13, 'W']],
        'Cuba': [[21, 'N'], [80, 'W']],
        'Cameroon': [[6, 'N'], [12, 'E']],
        'Montenegro': [[42, 'N'], [19, 'E']],
        'Republic of the Congo': [[1, 'S'], [15, 'E']],
        'Burkina Faso': [[13, 'N'], [2, 'W']],
        'Togo': [[8, 'N'], [1, 'E']],
        'Virgin Islands': [[18, 'N'], [64, 'W']],
        'China': [[35, 'N'], [105, 'E']],
        'Armenia': [[40, 'N'], [45, 'E']],
        'Timor-Leste': [[8, 'S'], [125, 'E']],
        'Dominican Republic': [[19, 'N'], [70, 'W']],
        'Ukraine': [[49, 'N'], [32, 'E']],
        'Bahrain': [[26, 'N'], [50, 'E']],
        'Tonga': [[20, 'S'], [175, 'W']],
        'Finland': [[64, 'N'], [26, 'E']],
        'Libya': [[25, 'N'], [17, 'E']],
        'Cayman Islands': [[19, 'N'], [80, 'W']],
        'Central African Republic': [[7, 'N'], [21, 'E']],
        'New Caledonia': [[21, 'S'], [165, 'E']],
        'Mauritius': [[20, 'S'], [57, 'E']],
        'Liechtenstein': [[47, 'N'], [9, 'E']],
        'Vietnam': [[16, 'N'], [107, 'E']],
        'British Virgin Islands': [[18, 'N'], [64, 'W']],
        'Mali': [[17, 'N'], [4, 'W']],
        'Vatican City': [[41, 'N'], [12, 'E']],
        'Russia': [[60, 'N'], [100, 'E']],
        'Bulgaria': [[43, 'N'], [25, 'E']],
        'United States': [[38, 'N'], [97, 'W']],
        'Romania': [[46, 'N'], [25, 'E']],
        'Angola': [[12, 'S'], [18, 'E']],
        'Chad': [[15, 'N'], [19, 'E']],
        'South Africa': [[29, 'S'], [24, 'E']],
        'Tokelau': [[9, 'S'], [172, 'W']],
        'Turks and Caicos Islands': [[21, 'N'], [71, 'W']],
        'South Georgia and the South Sandwich Islands': [[54, 'S'], [37, 'W']],
        'Sweden': [[62, 'N'], [15, 'E']],
        'Qatar': [[25, 'N'], [51, 'E']],
        'Malaysia': [[2, 'N'], [112, 'E']],
        'Senegal': [[14, 'N'], [14, 'W']],
        'Latvia': [[57, 'N'], [25, 'E']],
        'Clipperton Island': [[10, 'N'], [109, 'W']],
        'Uganda': [[1, 'N'], [32, 'E']],
        'Japan': [[36, 'N'], [138, 'E']],
        'Niger': [[16, 'N'], [8, 'E']],
        'Brazil': [[10, 'S'], [55, 'W']],
        'Faroe Islands': [[62, 'N'], [7, 'W']],
        'Guinea': [[11, 'N'], [10, 'W']],
        'Panama': [[9, 'N'], [80, 'W']],
        'Costa Rica': [[10, 'N'], [84, 'W']],
        'Luxembourg': [[49, 'N'], [6, 'E']],
        'Cape Verde': [[16, 'N'], [24, 'W']],
        'Andorra': [[42, 'N'], [1, 'E']],
        'Gibraltar': [[36, 'N'], [5, 'W']],
        'Ireland': [[53, 'N'], [8, 'W']],
        'Syria': [[35, 'N'], [38, 'E']],
        'Palau': [[7, 'N'], [134, 'E']],
        'Nigeria': [[10, 'N'], [8, 'E']],
        'Ecuador': [[2, 'S'], [77, 'W']],
        'Northern Mariana Islands': [[15, 'N'], [145, 'E']],
        'Brunei': [[4, 'N'], [114, 'E']],
        'Mozambique': [[18, 'S'], [35, 'E']],
        'Australia': [[27, 'S'], [133, 'E']],
        'Iran': [[32, 'N'], [53, 'E']],
        'Algeria': [[28, 'N'], [3, 'E']],
        'Svalbard': [[78, 'N'], [20, 'E']],
        'El Salvador': [[13, 'N'], [88, 'W']],
        'Tuvalu': [[8, 'S'], [178, 'E']],
        'Pitcairn Islands': [[25, 'S'], [130, 'W']],
        'Czech Republic': [[49, 'N'], [15, 'E']],
        'Marshall Islands': [[9, 'N'], [168, 'E']],
        'Chile': [[30, 'S'], [71, 'W']],
        'Puerto Rico': [[18, 'N'], [66, 'W']],
        'Belgium': [[50, 'N'], [4, 'E']],
        'Kiribati': [[1, 'N'], [173, 'E']],
        'Haiti': [[19, 'N'], [72, 'W']],
        'Belize': [[17, 'N'], [88, 'W']],
        'Hong Kong': [[22, 'N'], [114, 'E']],
        'Saint Lucia': [[13, 'N'], [60, 'W']],
        'Georgia': [[42, 'N'], [43, 'E']],
        'Mexico': [[23, 'N'], [102, 'W']],
        'Denmark': [[56, 'N'], [10, 'E']],
        'Poland': [[52, 'N'], [20, 'E']],
        'Moldova': [[47, 'N'], [29, 'E']],
        'Morocco': [[32, 'N'], [5, 'W']],
        'Namibia': [[22, 'S'], [17, 'E']],
        'Mongolia': [[46, 'N'], [105, 'E']],
        'Guernsey': [[49, 'N'], [2, 'W']],
        'Thailand': [[15, 'N'], [100, 'E']],
        'Switzerland': [[47, 'N'], [8, 'E']],
        'Grenada': [[12, 'N'], [61, 'W']],
        'Navassa Island': [[18, 'N'], [75, 'W']],
        'Isle of Man': [[54, 'N'], [4, 'W']],
        'Portugal': [[39, 'N'], [8, 'W']],
        'Estonia': [[59, 'N'], [26, 'E']],
        'Kosovo': [[42, 'N'], [21, 'E']],
        'Norfolk Island': [[29, 'S'], [167, 'E']],
        'Bouvet Island': [[54, 'S'], [3, 'E']],
        'Lebanon': [[33, 'N'], [35, 'E']],
        'Sierra Leone': [[8, 'N'], [11, 'W']],
        'Uzbekistan': [[41, 'N'], [64, 'E']],
        'Tunisia': [[34, 'N'], [9, 'E']],
        'Djibouti': [[11, 'N'], [43, 'E']],
        'Heard Island and McDonald Islands': [[53, 'S'], [72, 'E']],
        'Antigua and Barbuda': [[17, 'N'], [61, 'W']],
        'Spain': [[40, 'N'], [4, 'W']],
        'Colombia': [[4, 'N'], [72, 'W']],
        'Burundi': [[3, 'S'], [30, 'E']],
        'Taiwan': [[23, 'N'], [121, 'E']],
        'Cyprus': [[35, 'N'], [33, 'E']],
        'Barbados': [[13, 'N'], [59, 'W']],
        'Falkland Islands (Islas Malvinas)': [[51, 'S'], [59, 'W']],
        'Madagascar': [[20, 'S'], [47, 'E']],
        'Italy': [[42, 'N'], [12, 'E']],
        'Bhutan': [[27, 'N'], [90, 'E']],
        'Sudan': [[15, 'N'], [30, 'E']],
        'Vanuatu': [[16, 'S'], [167, 'E']],
        'Malta': [[35, 'N'], [14, 'E']],
        'Hungary': [[47, 'N'], [20, 'E']],
        'Democratic Republic of the Congo': [[0, 'N'], [25, 'E']],
        'Netherlands': [[52, 'N'], [5, 'E']],
        'Bermuda': [[32, 'N'], [64, 'W']],
        'Suriname': [[4, 'N'], [56, 'W']],
        'Anguilla': [[18, 'N'], [63, 'W']],
        'Venezuela': [[8, 'N'], [66, 'W']],
        'Netherlands Antilles': [[12, 'N'], [69, 'W']],
        'Israel': [[31, 'N'], [34, 'E']],
        'Paracel Islands': [[16, 'N'], [112, 'E']],
        'Wake Island': [[19, 'N'], [166, 'E']],
        'Indonesia': [[5, 'S'], [120, 'E']],
        'Iceland': [[65, 'N'], [18, 'W']],
        'Zambia': [[15, 'S'], [30, 'E']],
        'Samoa': [[13, 'S'], [172, 'W']],
        'Austria': [[47, 'N'], [13, 'E']],
        'Papua New Guinea': [[6, 'S'], [147, 'E']],
        'Malawi': [[13, 'S'], [34, 'E']],
        'Zimbabwe': [[20, 'S'], [30, 'E']],
        'Germany': [[51, 'N'], [9, 'E']],
        'Dhekelia': [[34, 'N'], [33, 'E']],
        'Kazakhstan': [[48, 'N'], [68, 'E']],
        'Philippines': [[13, 'N'], [122, 'E']],
        'Eritrea': [[15, 'N'], [39, 'E']],
        'Kyrgyzstan': [[41, 'N'], [75, 'E']],
        'Mayotte': [[12, 'S'], [45, 'E']],
        'Iraq': [[33, 'N'], [44, 'E']],
        'Montserrat': [[16, 'N'], [62, 'W']],
        'Coral Sea Islands': [[18, 'S'], [152, 'E']],
        'Macedonia': [[41, 'N'], [22, 'E']],
        'British Indian Ocean Territory': [[6, 'S'], [71, 'E']],
        'North Korea': [[40, 'N'], [127, 'E']],
        'Trinidad and Tobago': [[11, 'N'], [61, 'W']],
        'Akrotiri': [[34, 'N'], [32, 'E']],
        'Guyana': [[5, 'N'], [59, 'W']],
        'Belarus': [[53, 'N'], [28, 'E']],
        'Nepal': [[28, 'N'], [84, 'E']],
        'Burma': [[22, 'N'], [98, 'E']],
        'Honduras': [[15, 'N'], [86, 'W']],
        'Equatorial Guinea': [[2, 'N'], [10, 'E']],
        'Egypt': [[27, 'N'], [30, 'E']],
        'Nicaragua': [[13, 'N'], [85, 'W']],
        'Singapore': [[1, 'N'], [103, 'E']],
        'Serbia': [[44, 'N'], [21, 'E']],
        'Botswana': [[22, 'S'], [24, 'E']],
        'United Kingdom': [[54, 'N'], [2, 'W']],
        'Antarctica': [[90, 'S'], [0, 'E']],
        'Christmas Island': [[10, 'S'], [105, 'E']],
        'Greece': [[39, 'N'], [22, 'E']],
        'Sri Lanka': [[7, 'N'], [81, 'E']],
        'Croatia': [[45, 'N'], [15, 'E']],
        'Comoros': [[12, 'S'], [44, 'E']]
        }
    g = Graph()
    g.add_edges(edges)
    g.add_vertices(gps_coordinates)
    g.gps_coordinates = gps_coordinates
    g.name("World Map")
    return g
