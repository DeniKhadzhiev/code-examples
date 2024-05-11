# import sqlite3
# from flashcards import *
# import json
# import random

# def add_card(login : int, subject : str, term : str, meaning : str):
#     db_uq = sqlite3.connect('user_question.db')
#     cursor_db = db_uq.cursor()
#     cursor_db.execute('SELECT cards from flashcards WHERE login = ?', (login,))
#     req_res = cursor_db.fetchone()
#     d = json.loads(req_res)

#     if subject in d:
#         if term not in d[subject]:
#             d[subject][term] = meaning
#         else:
#             pass
#     else:
#         pass

#     cursor_db.execute('UPDATE flashcards SET cards = ? WHERE login = ?', (json.dumps(d), login))
#     db_uq.commit()
#     cursor_db.close()
#     db_uq.close()

    

# def add_subject(login : int, subject : str):
#     db_uq = sqlite3.connect('user_question.db')
#     cursor_db = db_uq.cursor()
#     cursor_db.execute('SELECT cards from flashcards WHERE login = ?', (login,))
#     req_res = cursor_db.fetchall()
#     d = json.loads(req_res)
#     if subject in d:
#         pass
#     else:
#         d[subject] = {}

#     cursor_db.execute('UPDATE flashcards SET cards = ? WHERE login = ?', (json.dumps(d), login))
#     db_uq.commit()
#     cursor_db.close()
#     db_uq.close()

# def delete_card(login : int, subject : str, term : str):
#     db_uq = sqlite3.connect('user_question.db')
#     cursor_db = db_uq.cursor()
#     cursor_db.execute('SELECT cards from flashcards WHERE login = ?', (login,))
#     req_res = cursor_db.fetchone()
#     d = json.loads(req_res)

#     if subject in d:
#         if term in d[subject]:
#             d[subject].pop(term)
#         else:
#             pass
#     else:
#         pass
    
#     cursor_db.execute('UPDATE flashcards SET cards = ? WHERE login = ?', (json.dumps(d), login))
#     db_uq.commit()
#     cursor_db.close()
#     db_uq.close()

# def delete_subject(login : int, subject : str):
#     db_uq = sqlite3.connect('user_question.db')
#     cursor_db = db_uq.cursor()
#     cursor_db.execute('SELECT cards from flashcards WHERE login = ?', (login,))
#     req_res = cursor_db.fetchone()
#     d = json.loads(req_res)

#     if subject in d:
#         d.pop(subject)
#     else:
#         pass
    
#     cursor_db.execute('UPDATE flashcards SET cards = ? WHERE login = ?', (json.dumps(d), login))
#     db_uq.commit()
#     cursor_db.close()
#     db_uq.close()



# def all_subjects(login : int) -> list:
#     db_uq = sqlite3.connect('user_question.db')
#     cursor_db = db_uq.cursor()
#     cursor_db.execute('SELECT cards from flashcards WHERE login = ?', (login,))
#     req_res = cursor_db.fetchall()
#     d = json.loads(req_res)
#     db_uq.commit()
#     cursor_db.close()
#     db_uq.close()

#     return list(d.keys())

# def all_subj_cards(login : int, subject : str) -> dict:
#     db_uq = sqlite3.connect('user_question.db')
#     cursor_db = db_uq.cursor()
#     cursor_db.execute('SELECT cards from flashcards WHERE login = ?', (login,))
#     req_res = cursor_db.fetchall()
#     d = json.loads(req_res)
#     db_uq.commit()
#     cursor_db.close()
#     db_uq.close()

#     if subject in list(d.keys()):
#         return d[subject] 
#     else:
#         pass

# def all_cards(login : int) -> dict:
#     db_uq = sqlite3.connect('user_question.db')
#     cursor_db = db_uq.cursor()
#     cursor_db.execute('SELECT cards from flashcards WHERE login = ?', (login,))
#     req_res = cursor_db.fetchall()
#     d = json.loads(req_res)
#     db_uq.commit()
#     cursor_db.close()
#     db_uq.close()

#     res = {}

#     for key in list(d.keys()):
#         res.update(d[key])
#     return res

# def flip_dict(d : dict) -> dict:
#     flipped = {j : i for i, j in d.items()}
#     return flipped

# def quiz(login: int, subject : str, quiz_mode : str, quiz_length : int) -> dict:   
#     # возможные значения quiz_mode : term, meaning
#     d = all_subj_cards(login, subject)
#     l = list(d.items())
#     random.shuffle(l)
#     d_shuffled = dict(l)
#     res = {}
#     m = min(quiz_length, len(l))
#     count = 0
#     if quiz_mode == "meaning":
#         for i in l:
#             if (count < m):
#                 res[i[0]] = i[1]
#             else:
#                 break
#     else:
#         for i in l:
#             if (count < m):
#                 res[i[1]] = i[0]
#             else:
#                 break
#     return {i : {"question" : key, "answer" : res[key]} for (i, key) in enumerate(res, start = 1)}
    

# # json.loads(input)

# def parse(inp : str) -> Flashcards:   
#     di = json.loads(inp)
#     for key in list(di.keys()):
#         # di[key] = list(map(lambda x: [Flashcard(x, j, x[j]) for j in x], di[key]))
#         for i in di[key]:
#             lst = []
#             tmp = Flashcard(key, i, di[key][i])
#             lst.append(tmp)
#         di[key] = lst
#     cards = Flashcards(di)
#     return cards


# db_uq = sqlite3.connect('user_question.db')
# cursor_db = db_uq.cursor()
# sql_create = '''CREATE TABLE flashcards(
# login TEXT PRIMARY KEY
# cards TEXT
# );'''

# #cursor_db.execute(sql_create)
# db_uq.commit()

# cursor_db.close()
# db_uq.close()



# # UPDATE db_uq SET cards = ... WHERE login = ...